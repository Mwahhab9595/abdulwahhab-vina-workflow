#!/usr/bin/env python3
"""Local server for Abdulwahhab Vina Workflow. Listens on 127.0.0.1 only.
Started by app.py (window) or by `python app.py --browser`."""
import csv, io, json, os, re, subprocess, sys, threading, time, zipfile
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs

FROZEN = bool(getattr(sys, "frozen", False))
HERE = getattr(sys, "_MEIPASS", None) or os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import version
from vina_backend import find_vina, vina_version
DATA = os.path.join(os.path.expanduser("~"), "AbdulwahhabVinaWorkflow")
RUNS = os.path.join(DATA, "runs")
SETTINGS_FILE = os.path.join(DATA, "settings.json")
PORT = int(os.environ.get("AVW_PORT", "8765"))
JOBS = {}
LOCK = threading.Lock()
DEFAULTS = {"exhaustiveness": 8, "num_poses": 10, "seed": 42, "cpu": 0, "scoring": "vina", "energy_range": 3.0,
            "min_rmsd": 1.0, "padding": 5.0, "min_box": 20.0, "max_box": 40.0, "rmsd_cutoff": 2.0, "ph": 7.4,
            "vina_path": ""}
LIMITS = {"exhaustiveness": (1, 128), "num_poses": (1, 20), "seed": (0, 2**31 - 1), "cpu": (0, 256),
          "energy_range": (0.5, 10), "min_rmsd": (0, 5), "padding": (0, 15), "min_box": (10, 60),
          "max_box": (20, 60), "rmsd_cutoff": (0.5, 5), "ph": (4, 10)}

def load_settings():
    s = dict(DEFAULTS)
    try: s.update(json.load(open(SETTINGS_FILE, encoding="utf-8")))
    except Exception: pass
    return s

def save_settings(p):
    s = load_settings()
    for k, d in DEFAULTS.items():
        if k not in p: continue
        v = p[k]
        if k in ("scoring",):
            if v not in ("vina", "vinardo"): raise ValueError("Scoring function must be vina or vinardo")
            s[k] = v
        elif k == "vina_path":
            v = str(v or "").strip()
            if v and not os.path.isfile(v): raise ValueError("The Vina program path you entered does not exist")
            s[k] = v
        else:
            try: v = type(d)(float(v)) if isinstance(d, float) else int(float(v))
            except Exception: raise ValueError(f"Setting '{k}' must be a number")
            lo, hi = LIMITS[k]
            if not lo <= v <= hi: raise ValueError(f"Setting '{k}' must be between {lo} and {hi}")
            s[k] = v
    os.makedirs(DATA, exist_ok=True)
    json.dump(s, open(SETTINGS_FILE, "w", encoding="utf-8"), indent=2)
    return s

def active_job():
    return next((j for j in JOBS.values() if j["state"] == "running"), None)

# ------------------------------------------------------------------ friendly errors
def friendly(code_text):
    t = code_text or ""
    rules = [
        (r"Vina program not found", "The Vina docking program was not found. Reinstall the app, or set its path in Settings."),
        (r"Vina failed", "Vina stopped with an error. The most common cause is a box that does not touch the protein or is too small. Check the box in step 2 and try again."),
        (r"protocol not validated", "Validation failed: the program could not put the known ligand back in its crystal position (RMSD above the limit). The results would not be trustworthy, so docking was stopped. Check the box and the chain choice, or tick 'Continue anyway (exploratory)' if you only want a rough look."),
        (r"box of .* exceeds|needs a box of", "The binding site is larger than the maximum box size. Choose a smaller site, select fewer residues, or raise the maximum box in Settings."),
        (r"no pocket found", "No binding pocket was found in this structure. Choose a structure that has a ligand, or draw the box by hand."),
        (r"None of the selected residues", "The residues you clicked are not in the selected chains. Select the chain that contains them."),
        (r"box centre/size must be", "The box numbers are not valid. Enter three numbers for the centre and three for the size, for example -1.8, 14.3, 16.7."),
        (r"Could not download", "The structure could not be downloaded from the PDB. Check the internet connection or upload the .pdb file."),
        (r"bad SMILES|embed failed", "One compound could not be read or converted to 3D. Check its SMILES or structure file."),
        (r"output directory", "The output folder is not empty. Start a new run."),
        (r"Traceback", "An unexpected error happened. Open the details below and send them to your instructor."),
    ]
    for pat, msg in rules:
        if re.search(pat, t): return msg
    return None

# ------------------------------------------------------------------ structure
PDB_CACHE = {}
def fetch_pdb_text(pdb, pdb_text):
    if (pdb_text or "").strip(): return pdb_text, "uploaded PDB file"
    pdb = (pdb or "").strip()
    if pdb.lower() == "example":
        return open(os.path.join(HERE, "examples", "3PTB.pdb"), encoding="utf-8").read(), "3PTB"
    if not re.fullmatch(r"[A-Za-z0-9]{4}", pdb):
        raise ValueError("Enter a 4-character PDB ID (for example 3PTB) or upload a .pdb file")
    pdb = pdb.upper()
    if pdb not in PDB_CACHE:
        import urllib.request
        try: PDB_CACHE[pdb] = urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb}.pdb", timeout=60).read().decode()
        except Exception as e:
            raise ValueError(f"Could not download {pdb} from the PDB ({e}). Check the PDB ID and the internet connection, or upload the .pdb file")
    return PDB_CACHE[pdb], pdb

def inspect_structure(pdb, pdb_text):
    import numpy as np
    from scipy.spatial import cKDTree
    import avw_analysis as A
    text, label = fetch_pdb_text(pdb, pdb_text)
    chains, atoms, hets = {}, [], {}
    for l in text.splitlines():
        if l.startswith("ENDMDL"): break
        if len(l) < 54: continue
        if l.startswith("ATOM") and l[16] in " A":
            c = l[21]; d = chains.setdefault(c, {"residues": set(), "atoms": 0}); d["residues"].add(l[22:27]); d["atoms"] += 1
            atoms.append((c, l[22:26].strip(), [float(l[30:38]), float(l[38:46]), float(l[46:54])]))
        elif l.startswith("HETATM") and l[16] in " A":
            r = l[17:20].strip()
            if r in A.NOT_LIGAND or l[76:78].strip() == "H": continue
            hets.setdefault((r, l[21], l[22:26].strip()), []).append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    if not chains: raise ValueError("No protein atoms were found in this structure. Is this a protein PDB file?")
    xyz = np.array([a[2] for a in atoms]); tree = cKDTree(xyz)
    ligs = []
    for (r, c, n), pts in hets.items():
        if len(pts) < 6: continue
        touch = {}
        for p in pts:
            for j in tree.query_ball_point(p, 4.0): touch.setdefault(atoms[j][0], set()).add(atoms[j][1])
        ligs.append({"res": r, "chain": c, "num": n, "heavy_atoms": len(pts), "contacts": {k: len(v) for k, v in sorted(touch.items())}})
    ligs.sort(key=lambda x: -x["heavy_atoms"])
    title = " ".join(l[10:].strip() for l in text.splitlines() if l.startswith("TITLE"))[:160]
    info = {"label": label, "title": title,
            "chains": [{"id": k, "residues": len(v["residues"]), "atoms": v["atoms"]} for k, v in sorted(chains.items())],
            "ligands": ligs[:12]}
    info["suggestion"] = A.suggest_chains(info, text)
    return info

def filter_chains(text, keep, het=True):
    out = []
    for l in text.splitlines():
        if l.startswith("ENDMDL"): break
        if len(l) > 21 and l[21] in keep:
            if l.startswith("ATOM"): out.append(l)
            elif het and l.startswith("HETATM") and l[17:20].strip() not in ("HOH", "WAT", "DOD"): out.append(l)
    return "\n".join(out) + "\nEND\n"

def prot_lines_for(p):
    text, _ = fetch_pdb_text(p.get("pdb"), p.get("pdb_text"))
    keep = set(p.get("chains") or [])
    if not keep: raise ValueError("Choose the protein chain(s) first (step 1)")
    lines = [l for l in filter_chains(text, keep, het=False).splitlines() if l.startswith("ATOM") and l[16] in " A"]
    if not lines: raise ValueError("The selected chains contain no protein atoms")
    return text, lines

def api_pockets(p):
    import dockflow as D, avw_analysis as A
    text, lines = prot_lines_for(p)
    return {"pockets": A.pocket_suggestions(lines, D.find_pockets, D.lining_residues)}

def api_box(p):
    import numpy as np, dockflow as D
    text, lines = prot_lines_for(p)
    st = load_settings(); pad = float(p.get("padding") or st["padding"])
    m = p.get("method"); note = None
    if m == "cocrystal":
        lig = (p.get("ligand_id") or "").strip(); keep = set(p.get("chains") or [])
        pts = []
        for l in text.splitlines():
            if l.startswith("HETATM") and l[17:20].strip() == lig and l[21] in keep and l[16] in " A" and l[76:78].strip() != "H":
                pts.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        if not pts: raise ValueError(f"Ligand {lig or '(none)'} was not found in the selected chain(s)")
        c, s, raw = D.box_from_points(np.array(pts), pad, st["min_box"])
    elif m == "residues":
        pts = D.residue_coords(lines, p.get("residues") or [])
        c, s, raw = D.box_from_points(pts, pad, st["min_box"])
    elif m == "pocket":
        import avw_analysis as A
        pk = D.find_pockets(np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines]))
        i = int(p.get("pocket_index") or 1)
        if i > len(pk): raise ValueError("That pocket does not exist")
        c, s, raw = D.box_from_points(pk[i - 1]["points"], pad, st["min_box"])
    elif m == "manual":
        try:
            c = np.array([float(x) for x in p["center"]]); s = np.array([float(x) for x in p["size"]]); raw = s
            assert len(c) == 3 and len(s) == 3 and (s > 0).all()
        except Exception: raise ValueError("Enter three numbers for the centre and three positive numbers for the size")
    else: raise ValueError("Unknown box method")
    if (raw > st["max_box"]).any() and m != "manual":
        note = f"This site needs a box of {', '.join(f'{x:.0f}' for x in raw)} A, above the maximum ({st['max_box']:.0f} A). Choose fewer residues or raise the limit in Settings."
    if m == "manual" and (s > 60).any(): raise ValueError("A manual box larger than 60 A is not allowed")
    return {"center": [round(float(x), 2) for x in c], "size": [round(float(x), 1) for x in s], "note": note}

# ------------------------------------------------------------------ jobs
def engine_cmd():
    if FROZEN: return [sys.executable, "--engine"]
    return [sys.executable, "-u", os.path.join(HERE, "app.py"), "--engine"]

def start_job(p):
    st = load_settings()
    comps = p.get("compounds") or []
    for c in comps:
        if c.get("level") == "reject": raise ValueError(f"'{c.get('name')}' was rejected by the input checks and cannot be docked")
        if c.get("level") == "warn" and not c.get("accepted"):
            raise ValueError(f"'{c.get('name')}' has input warnings that you have not accepted. Review its structure, then accept or remove it")
    rows = [(str(c.get("name", "")).strip(), str(c.get("smiles", "")).strip()) for c in comps if c.get("smiles")]
    rows = [r for r in rows if r[0] and r[1]]
    if not rows: raise ValueError("The compound list is empty. Add at least one compound in step 2 first")
    if len(rows) > 200: raise ValueError("Maximum 200 compounds per run")
    if not (st["vina_path"] or find_vina()): raise ValueError("The Vina docking program was not found. Reinstall the app or set its path in Settings")
    pdb = (p.get("pdb") or "").strip()
    text, label = fetch_pdb_text(pdb, p.get("pdb_text"))
    info = inspect_structure(pdb, p.get("pdb_text"))
    ids = [c["id"] for c in info["chains"]]
    keep = [c for c in (p.get("chains") or []) if c in ids]
    if not keep: raise ValueError("Choose which protein chain(s) to use in step 1")
    lig = (p.get("ligand_id") or "").strip()
    if lig:
        if not re.fullmatch(r"[A-Za-z0-9]{1,3}", lig): raise ValueError("Ligand ID is 1-3 letters/digits (for example BEN)")
        if not any(l["res"] == lig and l["chain"] in keep for l in info["ligands"]):
            raise ValueError(f"Ligand {lig} is not found in the selected chain(s) {','.join(keep)}. Pick the chain that carries it or change the Ligand ID")
    jid = time.strftime("%Y%m%d-%H%M%S"); base = os.path.join(RUNS, jid); os.makedirs(base)
    safe = re.sub(r"[^A-Za-z0-9]", "", label)[:8] or "input"
    pf = os.path.join(base, f"{safe}_chain{''.join(keep)}.pdb")
    open(pf, "w", encoding="utf-8").write(filter_chains(text, set(keep), het=True))
    lf = os.path.join(base, "library.csv")
    with open(lf, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["name", "smiles"])
        for r in rows: w.writerow(r)
    out = os.path.join(base, "results")
    with open(os.path.join(base, "input_check.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["name", "smiles", "source", "level", "warnings_accepted_by_user", "warnings", "note"])
        for c in comps:
            w.writerow([c.get("name"), c.get("smiles"), c.get("source", ""), c.get("level", ""), bool(c.get("accepted")),
                        " | ".join(c.get("warnings") or []), c.get("note", "")])
    def num(k, cast=float):
        v = p.get(k); return cast(float(v)) if v not in (None, "") else st[k]
    exh = num("exhaustiveness", int)
    if not 1 <= exh <= 128: raise ValueError("Exhaustiveness must be between 1 and 128")
    cmd = engine_cmd() + ["--pdb", pf, "--library", lf, "--out", out, "--exhaustiveness", str(exh),
          "--poses", str(num("num_poses", int)), "--seed", str(num("seed", int)), "--cpu", str(num("cpu", int)),
          "--scoring", p.get("scoring") or st["scoring"], "--energy-range", str(st["energy_range"]),
          "--min-rmsd", str(st["min_rmsd"]), "--padding", str(st["padding"]), "--min-box", str(st["min_box"]),
          "--max-box", str(st["max_box"]), "--rmsd-cutoff", str(st["rmsd_cutoff"]), "--ph", str(st["ph"])]
    if st["vina_path"]: cmd += ["--vina-path", st["vina_path"]]
    m = p.get("box_method") or "cocrystal"; boxdesc = ""
    if m == "cocrystal":
        if lig: cmd += ["--ligand-id", lig]
        boxdesc = "co-crystal ligand"
    elif m == "manual":
        c = p.get("box_center"); s = p.get("box_size")
        cmd += [f"--box-center={','.join(str(float(x)) for x in c)}", f"--box-size={','.join(str(float(x)) for x in s)}"]
        if lig: cmd += ["--ligand-id", lig]
        boxdesc = "manual coordinates"
    elif m == "residues":
        res = [r for r in (p.get("residues") or []) if re.fullmatch(r"[A-Za-z0-9]:-?\d+[A-Za-z]?", r)]
        if not res: raise ValueError("Click at least one residue in the 3D view to build the box")
        cmd += ["--box-residues", ",".join(res)]
        if lig: cmd += ["--ligand-id", lig]
        boxdesc = "selected residues"
    elif m == "pocket":
        cmd += ["--pocket-index", str(int(p.get("pocket_index") or 1))]
        if lig: cmd += ["--ligand-id", lig]
        boxdesc = "detected pocket"
    else: raise ValueError("Unknown box method")
    ref = (p.get("ref_smiles") or "").strip()
    if ref: cmd += ["--ref-smiles", ref, "--ref-name", (p.get("ref_name") or "reference").strip()[:40]]
    if p.get("allow_unvalidated"): cmd.append("--allow-unvalidated")
    job = {"id": jid, "state": "running", "log": [], "total": len(rows), "done": 0, "out": out,
           "target": f"{label} (chains {', '.join(keep)})", "started": time.time(), "code": None, "base": base,
           "selftest": bool(p.get("selftest"))}
    job["log"].append(f"[input] protein chains used (chosen by you): {', '.join(keep)} of {', '.join(ids)}; other chains are removed before preparation")
    job["log"].append(f"[input] box method: {boxdesc}; {len(rows)} compound(s); exhaustiveness {exh}")
    for l_ in info["ligands"]:
        if l_["chain"] in keep:
            out_c = [k for k in l_["contacts"] if k not in keep]
            if out_c: job["log"].append(f"[input] WARNING: ligand {l_['res']} also touches excluded chain(s) {', '.join(out_c)} within 4 A; part of its site is missing from this run")
    JOBS[jid] = job
    threading.Thread(target=run_proc, args=(job, cmd), daemon=True).start()
    return jid

def run_proc(job, cmd):
    try:
        kw = {"creationflags": 0x08000000} if os.name == "nt" else {}
        proc = subprocess.Popen(cmd, cwd=HERE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                                encoding="utf-8", errors="replace", **kw)
        job["proc"] = proc
        for line in proc.stdout:
            line = line.rstrip()
            if "[dockflow] docking " in line: job["done"] += 1
            line = line.replace(job["out"], "results").replace(os.path.expanduser("~"), "~")
            if line: job["log"].append(line)
        proc.wait(); job["code"] = proc.returncode
        if job["state"] != "cancelled": job["state"] = "finished" if proc.returncode == 0 else "failed"
    except Exception as e:
        job["log"].append(f"ERROR: {e}"); job["state"] = "failed"
    job["ended"] = time.time()

def status(jid):
    j = JOBS.get(jid)
    if not j: return None
    fm = friendly("\n".join(j["log"][-60:])) if j["state"] == "failed" else None
    return {"id": jid, "state": j["state"], "log": j["log"][-400:], "total": j["total"],
            "done": min(j["done"], j["total"]), "target": j["target"], "friendly": fm,
            "elapsed": int((j.get("ended") or time.time()) - j["started"])}

def rd(j, n):
    p = os.path.join(j["out"], n)
    return open(p, encoding="utf-8").read() if os.path.exists(p) else None

def results(jid):
    j = JOBS.get(jid)
    if not j: return None
    res = {"validation": None, "ranked": [], "report": "", "unvalidated": False, "selftest": j.get("selftest", False)}
    v = rd(j, "validation.json")
    if v: res["validation"] = json.loads(v)
    c = rd(j, "ranked.csv")
    if c: res["ranked"] = list(csv.DictReader(io.StringIO(c)))
    res["report"] = rd(j, "report.md") or ""
    res["unvalidated"] = os.path.exists(os.path.join(j["out"], "UNVALIDATED.txt"))
    return res

def diagram(jid, comp, pose):
    from rdkit import Chem
    import avw_analysis as A
    j = JOBS[jid]; store = json.loads(rd(j, "poses.json") or "{}")
    if comp not in store: raise ValueError("Compound not found in this run")
    ent = next((x for x in store[comp] if x["pose"] == pose), None)
    if not ent: raise ValueError("Pose not found")
    mol = Chem.MolFromMolBlock(ent["mol"])
    prot = [l for l in (rd(j, "receptor_used.pdb") or "").splitlines() if l.startswith("ATOM")]
    svg, inter = A.diagram_svg(mol, prot, f"{comp}, pose {pose}  (Vina {ent['vina']:.2f} kcal/mol)")
    return {"svg": svg, "interactions": inter}

def make_zip(jid):
    j = JOBS[jid]; buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for r, _, fs in os.walk(j["out"]):
            for f in fs:
                fp = os.path.join(r, f); z.write(fp, os.path.relpath(fp, j["out"]))
        ic = os.path.join(j["base"], "input_check.csv")
        if os.path.exists(ic): z.write(ic, "input_check.csv")
    return buf.getvalue()

def check_update():
    import urllib.request
    url = f"https://api.github.com/repos/{version.REPO}/releases/latest"
    try:
        d = json.loads(urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "AVW"}), timeout=15).read().decode())
    except Exception as e:
        return {"current": version.VERSION, "error": "Could not check for updates (no internet, or no release has been published yet)."}
    tag = (d.get("tag_name") or "").lstrip("vV")
    def tup(v): return tuple(int(x) if x.isdigit() else 0 for x in re.split(r"[.\-]", v)[:3])
    asset = next((a["browser_download_url"] for a in d.get("assets", []) if a["name"].lower().endswith((".zip", ".exe"))), None)
    return {"current": version.VERSION, "latest": tag, "newer": tup(tag) > tup(version.VERSION), "url": d.get("html_url"),
            "download": asset, "notes": (d.get("body") or "")[:1500]}

def about():
    import platform
    info = {"app": version.APP_NAME, "version": version.VERSION, "repo": version.REPO, "python": platform.python_version(),
            "vina": vina_version(load_settings()["vina_path"] or None) or "not found", "vina_path": load_settings()["vina_path"] or find_vina() or ""}
    try:
        from importlib.metadata import version as v
        info["libs"] = {k: v(k) for k in ("rdkit", "meeko", "openbabel-wheel", "numpy", "scipy") if True}
    except Exception: info["libs"] = {}
    f = os.path.join(HERE, "ABOUT.md")
    info["about_md"] = open(f, encoding="utf-8").read() if os.path.exists(f) else ""
    return info

class H(BaseHTTPRequestHandler):
    def log_message(self, *a): pass
    def _send(self, code, body, ctype="application/json; charset=utf-8", extra=None):
        if isinstance(body, (dict, list)): body = json.dumps(body, ensure_ascii=False)
        if isinstance(body, str): body = body.encode("utf-8")
        self.send_response(code); self.send_header("Content-Type", ctype); self.send_header("Content-Length", str(len(body)))
        for k, v in (extra or {}).items(): self.send_header(k, v)
        self.end_headers(); self.wfile.write(body)
    def _local_ok(self):
        return (self.headers.get("Host") or "").split(":")[0] in ("localhost", "127.0.0.1")
    def do_GET(self):
        if not self._local_ok(): return self._send(403, {"error": "forbidden"})
        u = urlparse(self.path); q = parse_qs(u.query); jid = (q.get("id") or [""])[0]
        try:
            if u.path == "/": return self._send(200, open(os.path.join(HERE, "gui", "index.html"), encoding="utf-8").read(), "text/html; charset=utf-8")
            if u.path.startswith("/vendor/"):
                f = os.path.join(HERE, "gui", "vendor", os.path.basename(u.path))
                if os.path.isfile(f): return self._send(200, open(f, "rb").read(), "application/javascript")
            if u.path == "/api/example":
                return self._send(200, {"library": open(os.path.join(HERE, "examples", "trypsin_library.csv"), encoding="utf-8").read()})
            if u.path == "/api/settings": return self._send(200, {"settings": load_settings(), "defaults": DEFAULTS, "limits": LIMITS})
            if u.path == "/api/about": return self._send(200, about())
            if u.path == "/api/update-check": return self._send(200, check_update())
            if u.path == "/api/status":
                s = status(jid); return self._send(200 if s else 404, s or {"error": "not found"})
            if u.path == "/api/results":
                r = results(jid); return self._send(200 if r else 404, r or {"error": "not found"})
            if u.path == "/api/poses" and jid in JOBS: return self._send(200, rd(JOBS[jid], "poses.json") or "{}")
            if u.path == "/api/receptor" and jid in JOBS: return self._send(200, rd(JOBS[jid], "receptor_used.pdb") or "", "text/plain")
            if u.path == "/api/diagram" and jid in JOBS:
                return self._send(200, diagram(jid, (q.get("compound") or [""])[0], int((q.get("pose") or ["1"])[0])))
            if u.path == "/api/zip" and jid in JOBS:
                return self._send(200, make_zip(jid), "application/zip", {"Content-Disposition": f'attachment; filename="AVW_{jid}.zip"'})
            if u.path == "/api/file" and jid in JOBS:
                n = (q.get("name") or [""])[0]
                if n in ("ranked.csv", "report.md", "best_poses.pdb", "all_poses_audit.csv", "validation.json"):
                    return self._send(200, rd(JOBS[jid], n) or "", "text/plain; charset=utf-8", {"Content-Disposition": f'attachment; filename="{n}"'})
            if u.path == "/api/jobs": return self._send(200, [{"id": k, "state": v["state"], "target": v["target"]} for k, v in JOBS.items()])
        except ValueError as e: return self._send(400, {"error": str(e)})
        except Exception as e: return self._send(500, {"error": f"Unexpected error: {e}"})
        self._send(404, {"error": "not found"})
    def do_POST(self):
        if not self._local_ok(): return self._send(403, {"error": "forbidden"})
        n = int(self.headers.get("Content-Length") or 0)
        if n > 40_000_000: return self._send(413, {"error": "That file is too large"})
        try: p = json.loads(self.rfile.read(n) or b"{}")
        except Exception: return self._send(400, {"error": "bad json"})
        u = urlparse(self.path)
        try:
            if u.path == "/api/run" or u.path == "/api/selftest":
                with LOCK:
                    if active_job(): return self._send(409, {"error": "A run is already in progress. Wait for it to finish or stop it."})
                    if u.path == "/api/selftest":
                        lib = list(csv.DictReader(open(os.path.join(HERE, "examples", "trypsin_library.csv"), encoding="utf-8")))[:1]
                        p = {"pdb": "example", "chains": ["A"], "ligand_id": "BEN", "box_method": "cocrystal", "selftest": True,
                             "exhaustiveness": 8, "compounds": [{"name": r["name"], "smiles": r["smiles"]} for r in lib]}
                    return self._send(200, {"id": start_job(p)})
            if u.path == "/api/parse":
                import ligand_input
                items = [(f.get("name"), f.get("text", "")) for f in (p.get("files") or [])]
                if (p.get("text") or "").strip(): items.insert(0, (None, p["text"]))
                return self._send(200, {"compounds": ligand_input.parse_many(items)})
            if u.path == "/api/inspect": return self._send(200, inspect_structure(p.get("pdb"), p.get("pdb_text")))
            if u.path == "/api/structure":
                text, _ = fetch_pdb_text(p.get("pdb"), p.get("pdb_text"))
                return self._send(200, {"pdb": filter_chains(text, set(p.get("chains") or []), het=True)})
            if u.path == "/api/pockets": return self._send(200, api_pockets(p))
            if u.path == "/api/box": return self._send(200, api_box(p))
            if u.path == "/api/settings": return self._send(200, {"settings": save_settings(p)})
            if u.path == "/api/cancel":
                j = JOBS.get(p.get("id"))
                if j and j["state"] == "running": j["state"] = "cancelled"; j["proc"].terminate()
                return self._send(200, {"ok": True})
        except ValueError as e: return self._send(400, {"error": str(e)})
        except Exception as e: return self._send(500, {"error": f"Unexpected error: {e}"})
        self._send(404, {"error": "not found"})

def make_server():
    os.makedirs(RUNS, exist_ok=True)
    for port in range(PORT, PORT + 20):
        try: return ThreadingHTTPServer(("127.0.0.1", port), H), port
        except OSError: continue
    raise RuntimeError("No free local port found")
