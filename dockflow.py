#!/usr/bin/env python3
"""
dockflow.py  -  automated AutoDock Vina pipeline with protocol validation
and multi-criteria pose review.

Steps
 1. Fetch/read target structure (PDB id or local file).
 2. Choose pocket: co-crystallised ligand if present, else --ref-smiles +
    geometric pocket finder (buriedness), else the top finder pocket.
 3. Build grid box from the pocket (ligand extent + padding).
 4. VALIDATE protocol: re-dock the co-crystallised ligand starting from a
    random conformer; PASS if the top pose RMSD <= --rmsd-cutoff (2.0 A).
 5. Dock the library (CSV: name,smiles).
 6. Review every pose on many criteria (not Vina score only) and pick the
    best pose per compound; write ranked table, report and PDB files.

Usage:
  python dockflow.py --pdb 1HSG --library lib.csv --out out
"""
import argparse, csv, json, os, sys, tempfile, urllib.request, math, time
from openbabel import pybel  # must be imported before rdkit/vina (swig clash)
import numpy as np
from scipy.spatial import cKDTree
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign, Descriptors
from rdkit.Chem import rdMolDescriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
from vina_backend import run_vina_cli

RDLogger.DisableLog("rdApp.*")
pybel.ob.obErrorLog.SetOutputLevel(0)

# residues that are NOT drug-like ligands
NOT_LIGAND = {"HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "MPD",
              "DMS", "ACT", "ACE", "NAG", "BMA", "MAN", "FMT", "TRS", "EPE", "MES",
              "CL", "NA", "K", "MG", "ZN", "CA", "MN", "FE", "CU", "NI", "CO", "CD",
              "IOD", "BR", "NO3", "SCN", "IMD", "BME", "CIT", "TAR", "MLI", "PG4", "1PE"}

DEFAULT_WEIGHTS = {  # sum is normalised automatically
    "vina": 0.30,        # Vina affinity
    "le": 0.10,          # ligand efficiency
    "interactions": 0.20,  # contacts with residues that bind the reference ligand
    "polar_contacts": 0.10,      # polar contacts with receptor
    "consensus": 0.10,   # reproducibility of the pose across the search
    "strain": 0.10,      # internal strain (MMFF), lower is better
    "buried": 0.10,      # fraction of ligand buried in the pocket
}

# ---------------------------------------------------------------- structure
def load_pdb(src, outdir):
    if os.path.isfile(src):
        return open(src).read(), os.path.basename(src).split(".")[0]
    url = f"https://files.rcsb.org/download/{src.upper()}.pdb"
    txt = urllib.request.urlopen(url, timeout=60).read().decode()
    open(os.path.join(outdir, f"{src.upper()}.pdb"), "w").write(txt)
    return txt, src.upper()

def parse(pdb_txt):
    prot, hets = [], {}
    for l in pdb_txt.splitlines():
        if l.startswith("ATOM"):
            if l[16] not in " A": continue          # keep altloc A only
            prot.append(l)
        elif l.startswith("HETATM"):
            if l[16] not in " A": continue
            res = l[17:20].strip()
            key = (res, l[21], l[22:26].strip())
            hets.setdefault(key, []).append(l)
    return prot, hets

def pick_ligand(hets, ligand_id=None):
    cands = []
    for (res, ch, num), lines in hets.items():
        if res in NOT_LIGAND: continue
        heavy = [x for x in lines if x[76:78].strip() != "H"]
        if ligand_id and res != ligand_id: continue
        if len(heavy) >= 6: cands.append(((res, ch, num), heavy))
    if not cands: return None
    cands.sort(key=lambda c: -len(c[1]))
    return cands[0]

def xyz(lines):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines])

def fetch_template(resname):
    try:
        url = f"https://files.rcsb.org/ligands/download/{resname}_ideal.sdf"
        txt = urllib.request.urlopen(url, timeout=60).read().decode()
        m = Chem.MolFromMolBlock(txt, removeHs=False)
        m = Chem.DeleteSubstructs(m, Chem.MolFromSmarts("[#1]"))
        Chem.SanitizeMol(m)
        return m
    except Exception:
        return None

# ---------------------------------------------------------------- pockets
def find_pockets(coords, spacing=1.5, ndir=30, top=3):
    tree = cKDTree(coords)
    lo, hi = coords.min(0) - 4, coords.max(0) + 4
    g = np.mgrid[lo[0]:hi[0]:spacing, lo[1]:hi[1]:spacing, lo[2]:hi[2]:spacing].reshape(3, -1).T
    d, _ = tree.query(g)
    g = g[(d > 3.2) & (d < 8)]                       # solvent-accessible shell
    i = np.arange(ndir) + 0.5
    phi = np.arccos(1 - 2 * i / ndir); th = np.pi * (1 + 5 ** 0.5) * i
    dirs = np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]
    hits = np.zeros(len(g), int)
    for dv in dirs:
        hit = np.zeros(len(g), bool)
        for step in np.arange(1.5, 11, 1.5):
            dd, _ = tree.query(g + dv * step)
            hit |= dd < 2.0
        hits += hit
    buried = g[hits >= int(ndir * 0.75)]
    if len(buried) == 0: return []
    bt = cKDTree(buried)
    lab = -np.ones(len(buried), int); n = 0
    for s in range(len(buried)):
        if lab[s] >= 0: continue
        stack = [s]; lab[s] = n
        while stack:
            for nb in bt.query_ball_point(buried[stack.pop()], spacing * 1.1):
                if lab[nb] < 0: lab[nb] = n; stack.append(nb)
        n += 1
    out = []
    for k in range(n):
        pts = buried[lab == k]
        out.append({"n_points": int(len(pts)), "center": pts.mean(0), "points": pts})
    out.sort(key=lambda p: -p["n_points"])
    return out[:top]


def residue_coords(prot_lines, residues):
    """residues: list of 'A:195' strings -> coordinates of all their atoms."""
    want = set()
    for r in residues:
        ch, num = r.split(":"); want.add((ch.strip(), num.strip()))
    sel = [l for l in prot_lines if (l[21], l[22:26].strip()) in want]
    if not sel: raise ValueError("None of the selected residues were found in the chosen chain(s)")
    return xyz(sel)

def lining_residues(prot_lines, pts, cutoff=4.5):
    tree = cKDTree(xyz(prot_lines)); idx = set()
    for p in pts: idx.update(tree.query_ball_point(p, cutoff))
    seen = {}
    for i in sorted(idx):
        l = prot_lines[i]; seen[(l[21], l[22:26].strip())] = l[17:20] + l[22:26].strip() + ":" + l[21]
    return [seen[k] for k in sorted(seen, key=lambda k: (k[0], int(k[1]) if k[1].lstrip("-").isdigit() else 0))]

def box_from_points(pts, padding, min_box, max_box=None):
    ext = pts.max(0) - pts.min(0); center = (pts.max(0) + pts.min(0)) / 2
    raw = ext + 2 * padding
    return center, np.maximum(raw, min_box), raw

# ---------------------------------------------------------------- receptor
def prep_receptor(prot_lines, outdir, ph=7.4):
    pdb = os.path.join(outdir, "receptor_clean.pdb")
    open(pdb, "w").write("\n".join(prot_lines) + "\nEND\n")
    mol = pybel.readstring("pdb", open(pdb).read())
    mol.OBMol.AddHydrogens(False, True, ph)
    qt = os.path.join(outdir, "receptor.pdbqt")
    mol.write("pdbqt", qt, overwrite=True, opt={"r": None})
    return pdb, qt

# ---------------------------------------------------------------- ligand
def mol_from_smiles(smi, seed=7):
    m = Chem.MolFromSmiles(smi)
    if m is None: raise ValueError(f"bad SMILES: {smi}")
    m = Chem.AddHs(m)
    ps = AllChem.ETKDGv3(); ps.randomSeed = seed; ps.useRandomCoords = True
    if AllChem.EmbedMolecule(m, ps) != 0:
        raise ValueError(f"embed failed: {smi}")
    try: AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except Exception: pass
    return m

def to_pdbqt(m):
    prep = MoleculePreparation()
    setups = prep.prepare(m)
    s, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok: raise ValueError(err)
    return s

VINA_OPTS = {"seed": 42, "cpu": 0, "scoring": "vina", "energy_range": 3.0, "min_rmsd": 1.0, "exe": None}

def run_vina(rec_qt, lig_qt, center, size, exh, npose=10, seed=None):
    o = dict(VINA_OPTS)
    if seed is not None: o["seed"] = seed
    txt, en = run_vina_cli(rec_qt, lig_qt, [float(x) for x in center], [float(x) for x in size],
                           exh=exh, npose=npose, **o)
    pm = PDBQTMolecule(txt, is_dlg=False, skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pm, only_cluster_leads=False, keep_flexres=False)
    return mols[0], en

def heavy(m):
    return Chem.RemoveHs(m)

def conf_mol(m, cid):
    c = Chem.Mol(m, confId=cid) if False else Chem.Mol(m)
    c.RemoveAllConformers()
    c.AddConformer(Chem.Conformer(m.GetConformer(cid)), assignId=True)
    return c

def rmsd_sym(a, b):
    """symmetry-aware RMSD, no re-alignment."""
    try:
        return float(rdMolAlign.CalcRMS(a, b))
    except Exception:
        return float("nan")

# ---------------------------------------------------------------- metrics
class Receptor:
    def __init__(self, prot_lines, key_residues):
        self.lines = prot_lines
        self.xyz = xyz(prot_lines)
        self.el = np.array([l[76:78].strip() or l[12:16].strip()[0] for l in prot_lines])
        self.res = [(l[21], l[22:26].strip(), l[17:20]) for l in prot_lines]
        self.tree = cKDTree(self.xyz)
        self.polar = np.isin(self.el, ["N", "O"])
        self.key = key_residues
        self.res_idx = {}
        for i, r in enumerate(self.res):
            self.res_idx.setdefault(r[:2], []).append(i)

    def key_from_ligand(self, lig_xyz, cutoff=4.0):
        idx = set()
        for p in lig_xyz:
            idx.update(self.tree.query_ball_point(p, cutoff))
        return {(self.res[i][0], self.res[i][1]) for i in idx}

    def pose_metrics(self, pm):
        """pm: heavy-atom rdkit mol with 1 conformer."""
        pos = pm.GetConformer().GetPositions()
        el = [a.GetSymbol() for a in pm.GetAtoms()]
        n = len(pos)
        # clashes: any pair < 2.2 A, or non-polar pair < 3.0 A
        clash = 0
        for i, p in enumerate(pos):
            for j in self.tree.query_ball_point(p, 3.0):
                d = np.linalg.norm(self.xyz[j] - p)
                polar_pair = el[i] in ("N", "O") and self.polar[j]
                if d < 2.2 or (not polar_pair and d < 3.0 and d < 3.0 - 0.0 and d < 2.9):
                    clash += 1
        # hbond-like polar contacts
        hb = 0
        for i, p in enumerate(pos):
            if el[i] in ("N", "O"):
                for j in self.tree.query_ball_point(p, 3.4):
                    if self.polar[j] and np.linalg.norm(self.xyz[j] - p) >= 2.4:
                        hb += 1
        # burial: ligand heavy atoms with >= 10 receptor atoms within 4.5 A
        bur = sum(len(self.tree.query_ball_point(p, 4.5)) >= 10 for p in pos) / n
        # key residue recovery
        touched = set()
        for p in pos:
            for j in self.tree.query_ball_point(p, 4.0):
                touched.add(self.res[j][:2])
        key = (len(touched & {k[:2] for k in self.key}) / len(self.key)) if self.key else None
        return {"clashes": clash, "polar_contacts": hb, "buried": bur, "key_recovery": key}

def min_conf_energy(mol_noH, nconf=40):
    """lowest MMFF94 energy among random conformers (reference for strain)."""
    try:
        ref = Chem.AddHs(Chem.RemoveHs(mol_noH))
        if not AllChem.MMFFHasAllMoleculeParams(ref): return None
        ps = AllChem.ETKDGv3(); ps.randomSeed = 11
        AllChem.EmbedMultipleConfs(ref, nconf, ps)
        res = AllChem.MMFFOptimizeMoleculeConfs(ref, maxIters=1000)
        return min(e for _, e in res)
    except Exception:
        return None

def strain_energy(pose, e_min):
    """MMFF94 energy of the locally relaxed pose minus the lowest-energy conformer found (kcal/mol, rough)."""
    try:
        if e_min is None: return None
        m = Chem.AddHs(Chem.RemoveHs(pose), addCoords=True)
        props = AllChem.MMFFGetMoleculeProperties(m)
        ff = AllChem.MMFFGetMoleculeForceField(m, props)
        ff.Minimize(maxIts=500)
        return max(0.0, ff.CalcEnergy() - e_min)
    except Exception:
        return None

# ---------------------------------------------------------------- scoring
def clip01(x): return max(0.0, min(1.0, x))

def composite(row, weights, has_key):
    w = dict(weights)
    if not has_key: w["interactions"] = 0
    comp = {
        "vina": clip01((-row["vina"] - 4.0) / 8.0),            # -4..-12 kcal/mol
        "le": clip01((row["le"] - 0.15) / 0.30),                # 0.15..0.45
        "interactions": row["key_recovery"] if row["key_recovery"] is not None else 0,
        "polar_contacts": clip01(row["polar_contacts"] / 4.0),
        "consensus": row["consensus"],
        "strain": 1.0 - clip01((row["strain"] if row["strain"] is not None else 10) / 20.0),
        "buried": row["buried"],
    }
    tot = sum(w.values())
    s = sum(w[k] * comp[k] for k in w) / tot
    if row["clashes"] > 0: s *= 0.5                              # steric clash penalty
    return s, comp

def pose_pdb_block(m, resname="LIG", name="lig"):
    mh = Chem.AddHs(Chem.RemoveHs(m), addCoords=True)
    b = Chem.MolToPDBBlock(mh)
    out = []
    for l in b.splitlines():
        if l.startswith(("HETATM", "ATOM")):
            l = "HETATM" + l[6:17] + resname[:3].ljust(3) + l[20:]
            out.append(l)
    return "\n".join(out)

# ---------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", required=True, help="PDB id or local .pdb file")
    ap.add_argument("--library", required=True, help="CSV with columns name,smiles")
    ap.add_argument("--out", default="dock_out")
    ap.add_argument("--ligand-id", help="HET code of co-crystal ligand (default: largest)")
    ap.add_argument("--ref-smiles", help="literature reference ligand (when no co-crystal ligand)")
    ap.add_argument("--ref-name", default="reference")
    ap.add_argument("--padding", type=float, default=5.0, help="A added each side of the pocket ligand")
    ap.add_argument("--min-box", type=float, default=20.0)
    ap.add_argument("--max-box", type=float, default=40.0)
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--poses", type=int, default=10)
    ap.add_argument("--rmsd-cutoff", type=float, default=2.0)
    ap.add_argument("--top", type=int, default=5, help="compounds written to best_poses.pdb")
    ap.add_argument("--box-center", help="manual box centre 'x,y,z' (box method: manual coordinates)")
    ap.add_argument("--box-size", help="manual box size 'x,y,z' in A (used with --box-center)")
    ap.add_argument("--box-residues", help="box around these residues, e.g. 'A:195,A:189,A:57' (box method: residues)")
    ap.add_argument("--pocket-index", type=int, help="use pocket N (1-3) from the pocket finder (box method: pocket detection)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--cpu", type=int, default=0, help="0 = all cores")
    ap.add_argument("--scoring", default="vina", choices=["vina", "vinardo"])
    ap.add_argument("--energy-range", type=float, default=3.0)
    ap.add_argument("--min-rmsd", type=float, default=1.0)
    ap.add_argument("--ph", type=float, default=7.4, help="pH for receptor protonation")
    ap.add_argument("--vina-path", help="path to the Vina program (default: bundled copy)")
    ap.add_argument("--weights", help="JSON file overriding score weights")
    ap.add_argument("--retry-validation", action="store_true", help="if redock fails, retry once at 4x exhaustiveness (slow)")
    ap.add_argument("--allow-unvalidated", action="store_true",
                    help="EXPLORATORY: continue even if the protocol is not validated (outputs are labelled UNVALIDATED)")
    a = ap.parse_args(argv)
    VINA_OPTS.update(seed=a.seed, cpu=a.cpu, scoring=a.scoring, energy_range=a.energy_range, min_rmsd=a.min_rmsd, exe=a.vina_path)

    if os.path.isdir(a.out) and os.listdir(a.out):
        sys.exit(f"ERROR: output directory '{a.out}' is not empty. Use a new directory so stale "
                 "ranked/best-pose files from an earlier run can never be mistaken for this run.")
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    log = lambda *s: print("[dockflow]", *s, flush=True)
    weights = DEFAULT_WEIGHTS.copy()
    if a.weights: weights.update(json.load(open(a.weights)))

    pdb_txt, pid = load_pdb(a.pdb, a.out)
    prot, hets = parse(pdb_txt)
    rec_pdb, rec_qt = prep_receptor(prot, a.out, a.ph)
    log(f"target {pid}: {len(prot)} protein atoms")

    val = {"target": pid, "protocol_validated": False}
    xtal_mol = None; pocket_lig_xyz = None; ref_mol = None
    pick = pick_ligand(hets, a.ligand_id)
    if pick:
        (res, ch, num), llines = pick
        pocket_lig_xyz = xyz(llines)
        val.update(pocket_source=f"co-crystal ligand {res} chain {ch} {num}")
        tmpl = fetch_template(res)
        blk = "\n".join(llines) + "\nEND\n"
        raw = Chem.MolFromPDBBlock(blk, removeHs=True, proximityBonding=True)
        if tmpl is not None and raw is not None:
            for variant in ("as-is", "neutralised"):
                try:
                    t = Chem.Mol(tmpl)
                    if variant == "neutralised":
                        for at in t.GetAtoms():
                            at.SetFormalCharge(0); at.SetNoImplicit(False); at.SetNumExplicitHs(0)
                        Chem.SanitizeMol(t)
                    xtal_mol = AllChem.AssignBondOrdersFromTemplate(t, raw)
                    break
                except Exception as e:
                    log(f"template bond-order assignment ({variant}) failed:", e)
        ref_smiles = Chem.MolToSmiles(xtal_mol) if xtal_mol is not None else None
    elif a.ref_smiles:
        ref_smiles = a.ref_smiles
    else:
        ref_smiles = None

    # ---- pocket / box
    box_note = None
    if a.box_center:
        try:
            center = np.array([float(x) for x in a.box_center.split(",")])
            size = np.array([float(x) for x in (a.box_size or "20,20,20").split(",")])
            assert len(center) == 3 and len(size) == 3 and (size > 0).all()
        except Exception:
            sys.exit("ERROR: box centre/size must be three numbers each, like 10.5,-3,22")
        if (size > 60).any(): sys.exit("ERROR: manual box larger than 60 A is not allowed (slow and unreliable).")
        pocket_lig_xyz = None; val["pocket_source"] = "manual box coordinates"
        key_pts = None
    elif a.box_residues:
        pts = residue_coords(prot, [x for x in a.box_residues.split(",") if x.strip()])
        center, size, raw = box_from_points(pts, a.padding, a.min_box)
        if (raw > a.max_box).any():
            sys.exit(f"ERROR: the selected residues need a box of {np.round(raw,1)} A, above the maximum {a.max_box}. Select fewer residues or raise the limit in Settings.")
        pocket_lig_xyz = pts; key_pts = pts
        val["pocket_source"] = "box around selected residues: " + a.box_residues
    elif a.pocket_index:
        log("pocket finder")
        pk = find_pockets(xyz(prot))
        if not pk or a.pocket_index > len(pk): sys.exit("ERROR: pocket %d not found" % a.pocket_index)
        val["pockets_found"] = [{"rank": i + 1, "n_points": p["n_points"],
                                  "center": [round(float(x), 2) for x in p["center"]]} for i, p in enumerate(pk)]
        pts = pk[a.pocket_index - 1]["points"]
        center, size, raw = box_from_points(pts, a.padding, a.min_box)
        if (raw > a.max_box).any():
            size = np.minimum(raw, a.max_box); box_note = "pocket larger than the maximum box; box limited to %.0f A" % a.max_box
        pocket_lig_xyz = pts; key_pts = pts
        val["pocket_source"] = f"pocket finder, pocket {a.pocket_index} (geometric; verify visually)"
    else:
        key_pts = None
        if pocket_lig_xyz is None:
            log("no co-crystal ligand: running pocket finder")
            pk = find_pockets(xyz(prot))
            if not pk:
                sys.exit("no pocket found; supply a ligand-bound structure")
            val["pockets_found"] = [{"rank": i + 1, "n_points": p["n_points"],
                                      "center": [round(float(x), 2) for x in p["center"]]} for i, p in enumerate(pk)]
            pocket_lig_xyz = pk[0]["points"]
            val.update(pocket_source="geometric pocket finder (largest buried cavity); verify manually")
        center, size, raw_size = box_from_points(pocket_lig_xyz, a.padding, a.min_box)
        if (raw_size > a.max_box).any():
            sys.exit(f"ERROR: pocket/reference needs a box of {np.round(raw_size,1)} A which exceeds --max-box "
                     f"{a.max_box}. The box is NOT silently clipped (it would truncate the reference). "
                     "Raise --max-box, lower --padding, or choose another pocket.")
        assert np.all(pocket_lig_xyz >= center - size / 2) and np.all(pocket_lig_xyz <= center + size / 2)
    val["grid_center"] = [round(float(x), 3) for x in center]
    val["grid_size"] = [round(float(x), 2) for x in size]
    if box_note: val["box_note"] = box_note
    log("box center", val["grid_center"], "size", val["grid_size"])

    rec = Receptor(prot, set())
    key = rec.key_from_ligand(pocket_lig_xyz, 4.0) if pocket_lig_xyz is not None else set()
    rec.key = key
    val["key_residues"] = sorted(f"{r[2] if len(r)>2 else ''}{r[0]}:{r[1]}" for r in key)
    # readable key residue names
    names = {}
    for l in prot:
        names[(l[21], l[22:26].strip())] = l[17:20] + l[22:26].strip() + ":" + l[21]
    val["key_residues"] = sorted(names[k] for k in key if k in names)

    # ---- validation by re-docking
    if xtal_mol is not None and not (np.all(xyz(llines) >= center - size / 2 - 0.5) and np.all(xyz(llines) <= center + size / 2 + 0.5)):
        val["validation_note"] = ("The co-crystal ligand is outside the chosen box, so redocking validation was not possible with this box. "
                                "Use the co-crystal box method to validate the protocol.")
        xtal_mol = None
    if xtal_mol is not None:
        log("validating: re-docking", ref_smiles)
        try:
            start = mol_from_smiles(Chem.MolToSmiles(xtal_mol), seed=3)
            hx = heavy(xtal_mol)
            for exh in ((a.exhaustiveness, a.exhaustiveness * 4) if a.retry_validation else (a.exhaustiveness,)):
                m, en = run_vina(rec_qt, to_pdbqt(start), center, size, exh, a.poses)
                rm = []
                for cid in range(m.GetNumConformers()):
                    rm.append(rmsd_sym(heavy(conf_mol(m, cid)), hx))
                val.update(redock_exhaustiveness=exh, redock_top_rmsd=round(rm[0], 2),
                           redock_best_rmsd_any_pose=round(float(np.nanmin(rm)), 2),
                           redock_top_score=round(en[0], 2),
                           redock_rmsds=[round(x, 2) for x in rm],
                           redock_scores=[round(x, 2) for x in en])
                log(f"redock exh={exh}: top RMSD {rm[0]:.2f}  best-any {np.nanmin(rm):.2f}")
                if rm[0] <= a.rmsd_cutoff: break
            val["protocol_validated"] = bool(rm[0] <= a.rmsd_cutoff)
            val["viewer"] = {"crystal_pdb": pose_pdb_block(xtal_mol, 'XTL'), "redock_pdb": pose_pdb_block(conf_mol(m, 0), 'RDK')}
            val["validation_rule"] = f"top-ranked pose RMSD <= {a.rmsd_cutoff} A to crystal pose"
            ref_mol = xtal_mol
        except Exception as e:
            val["validation_error"] = str(e)
            log("validation error:", e)
    elif "validation_note" not in val:
        val["validation_note"] = ("No co-crystal pose available: RMSD validation impossible. "
                                  "Reference ligand (if given) is docked as a positive control only.")
    val["status"] = "VALIDATED" if val["protocol_validated"] else "UNVALIDATED"
    if not val["protocol_validated"]:
        val["gate"] = ("Re-docking validation failed or was impossible (no experimental pose). "
                       + ("Continuing only because --allow-unvalidated was given." if a.allow_unvalidated
                          else "Stopped: no screening, ranking or best-pose export was produced."))
        json.dump(val, open(os.path.join(a.out, "validation.json"), "w"), indent=2)
        if not a.allow_unvalidated:
            sys.exit("STOP: protocol not validated (see validation.json). Use --allow-unvalidated only for exploratory runs.")
        open(os.path.join(a.out, "UNVALIDATED.txt"), "w").write(
            "Outputs in this folder come from a protocol that was NOT validated against an experimental pose. Exploratory only.\n")

    # ---- library
    lib = list(csv.DictReader(open(a.library)))
    if a.ref_smiles and not pick:
        lib.insert(0, {"name": a.ref_name + " (control)", "smiles": a.ref_smiles})
    rows, best_blocks, best_mols, all_pose_rows = [], {}, {}, []
    pose_store = {}
    for it in lib:
        nm, smi = it["name"], it["smiles"]
        log("docking", nm)
        try:
            lm = mol_from_smiles(smi)
            m, en = run_vina(rec_qt, to_pdbqt(lm), center, size, a.exhaustiveness, a.poses)
        except Exception as e:
            rows.append({"compound": nm, "smiles": smi, "error": str(e)}); log("  failed:", e); continue
        nh = lm.GetNumHeavyAtoms()
        poses = [heavy(conf_mol(m, c)) for c in range(m.GetNumConformers())]
        best_e = en[0]
        e_min = min_conf_energy(poses[0])
        cand = []
        near = [j for j in range(len(poses)) if en[j] <= best_e + 2.0]
        for k, (pm, e) in enumerate(zip(poses, en)):
            # consensus: share of low-energy poses (within 2 kcal) that sit within 2 A of this pose
            cons = np.mean([rmsd_sym(pm, poses[j]) <= 2.0 for j in near])
            met = rec.pose_metrics(pm)
            r = {"compound": nm, "smiles": smi, "pose": k + 1, "vina": e, "le": -e / nh,
                 "consensus": float(cons), "heavy_atoms": nh, **met}
            r["strain"] = strain_energy(pm, e_min)
            r["score"], r["components"] = composite(r, weights, bool(rec.key))
            cand.append((r, pm))
        all_pose_rows.extend(r for r, _ in cand)       # audit: EVERY pose
        pose_store[nm] = [{"pose": r["pose"], "vina": round(r["vina"], 3), "score": round(r["score"], 3),
                           "polar_contacts": r["polar_contacts"], "clashes": r["clashes"], "pdb": pose_pdb_block(pm),
                           "mol": Chem.MolToMolBlock(pm)}
                          for r, pm in sorted(cand, key=lambda c: c[0]["pose"])]
        cand.sort(key=lambda c: -c[0]["score"])
        best, bm = cand[0]
        best["vina_top_pose"] = en[0]
        best["best_by_vina_is_best_overall"] = (best["pose"] == 1)
        rows.append(best); best_mols[nm] = bm
        log(f"  best pose #{best['pose']} vina {best['vina']:.2f} composite {best['score']:.3f}")

    ok = [r for r in rows if "error" not in r]
    ok.sort(key=lambda r: -r["score"])
    for i, r in enumerate(ok): r["rank"] = i + 1

    # ---- outputs
    cols = ["rank", "compound", "pose", "score", "vina", "vina_top_pose", "le", "key_recovery",
            "polar_contacts", "consensus", "strain", "buried", "clashes", "heavy_atoms", "smiles"]
    with open(os.path.join(a.out, "ranked.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(cols)
        for r in ok:
            w.writerow([round(r[c], 3) if isinstance(r.get(c), float) else r.get(c) for c in cols])
    acols = ["compound", "pose", "vina", "le", "key_recovery", "polar_contacts", "consensus", "strain",
             "buried", "clashes", "score"]
    with open(os.path.join(a.out, "all_poses_audit.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(acols + ["chosen_best_for_compound"])
        chosen = {(r["compound"], r["pose"]) for r in ok}
        for r in all_pose_rows:
            w.writerow([round(r[c], 3) if isinstance(r.get(c), float) else r.get(c) for c in acols]
                       + [int((r["compound"], r["pose"]) in chosen)])
    json.dump(ok, open(os.path.join(a.out, "ranked.json"), "w"), indent=2, default=float)
    json.dump(pose_store, open(os.path.join(a.out, "poses.json"), "w"))
    open(os.path.join(a.out, "receptor_used.pdb"), "w").write("\n".join(prot) + "\nEND\n")

    cdir = os.path.join(a.out, "complexes"); os.makedirs(cdir, exist_ok=True)
    rec_block = "\n".join(l for l in prot)
    with open(os.path.join(a.out, "best_poses.pdb"), "w") as f:
        f.write(f"REMARK dockflow {val['status']} protocol; heuristic shortlist for {pid}; one MODEL per compound\n")
        for r in ok[:a.top]:
            f.write(f"MODEL {r['rank']}\n")
            f.write(f"REMARK {r['compound']} pose {r['pose']} vina {r['vina']:.2f} composite {r['score']:.3f}\n")
            f.write(pose_pdb_block(best_mols[r["compound"]]) + "\nENDMDL\n")
        f.write("END\n")
    for r in ok:
        safe = "".join(c if c.isalnum() else "_" for c in r["compound"])
        with open(os.path.join(cdir, f"{r['rank']:02d}_{safe}_complex.pdb"), "w") as f:
            f.write(f"REMARK {r['compound']} best pose {r['pose']} vina {r['vina']:.2f} composite {r['score']:.3f}\n")
            f.write(rec_block + "\nTER\n" + pose_pdb_block(best_mols[r["compound"]]) + "\nEND\n")

    val["elapsed_s"] = round(time.time() - t0)
    val["weights"] = weights
    json.dump(val, open(os.path.join(a.out, "validation.json"), "w"), indent=2)

    # ---- report
    L = [f"# Docking report: {pid}  [{val['status']} protocol]", ""] + (["**EXPLORATORY: protocol not validated against an experimental pose.**", ""] if not val["protocol_validated"] else []) + [
         "## 1. Protocol validation", f"- Pocket: {val.get('pocket_source')}",
         f"- Grid center {val['grid_center']}, size {val['grid_size']} A"]
    if "redock_top_rmsd" in val:
        L += [f"- Re-docking of co-crystal ligand: top pose RMSD **{val['redock_top_rmsd']} A**, "
              f"best of {len(val['redock_rmsds'])} poses {val['redock_best_rmsd_any_pose']} A, score {val['redock_top_score']}",
              f"- Result: **{'PASS' if val['protocol_validated'] else 'FAIL - treat results with caution'}** "
              f"({val['validation_rule']})"]
    else:
        L.append(f"- {val.get('validation_note')}")
    L += ["", "## 2. Heuristic shortlist (multi-criteria; weights are unbenchmarked defaults)", "",
          "| Rank | Compound | Pose | Composite | Vina | LE | Key-res | Polar contacts | Consensus | Strain | Clash |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    f2 = lambda x: "-" if x is None else f"{x:.2f}"
    for r in ok:
        L.append(f"| {r['rank']} | {r['compound']} | {r['pose']} | {r['score']:.3f} | {r['vina']:.2f} | {r['le']:.2f} | "
                 f"{f2(r['key_recovery'])} | {r['polar_contacts']} | {r['consensus']:.2f} | {f2(r['strain'])} | {r['clashes']} |")
    diff = [r for r in ok if not r["best_by_vina_is_best_overall"]]
    L += ["", f"- Compounds where the chosen pose is NOT the top-Vina pose: {len(diff)}"
          + (" (" + ", ".join(r['compound'] for r in diff) + ")" if diff else ""),
          "", "Composite = weighted mix of: " + ", ".join(f"{k} {v}" for k, v in weights.items())
          + "; x0.5 if steric clash. Weights are adjustable (--weights).",
          "", "## 3. Limits",
          "- Vina scores are estimates, not binding free energies; do not compare across targets.",
          "- Receptor rigid; ligands exactly as written in SMILES (no pH enumeration).",
          "- Strain = MMFF94 vs best of 40 random conformers; a rough estimate.",
          "- Polar contacts = N/O...N/O distance 2.4-3.4 A only; not proof of hydrogen bonds (no angle/donor-acceptor check).",
          "- Every pose is scored in all_poses_audit.csv; the shortlist is a heuristic, not a proof of the best pose.",
          "- Receptor prep is automatic but approximate: protonation/tautomers, biological assembly, essential waters/metals/cofactors, ambiguous ligand copies and the choice of site need human review.",
          "- Wet-lab/MD confirmation is needed before conclusions."]
    open(os.path.join(a.out, "report.md"), "w").write("\n".join(L) + "\n")
    log("done in", val["elapsed_s"], "s ->", a.out)

if __name__ == "__main__":
    main()
