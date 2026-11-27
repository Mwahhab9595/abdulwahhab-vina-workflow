"""Ligand input reader for the GUI.

Turns SMILES text/CSV/.smi, MOL, SDF (multi-molecule) and MOL2 (multi-molecule)
into a list of (name, canonical SMILES) for dockflow.py, which only reads SMILES.
What is kept: atoms, bonds, formal charges and stereo that the file defines.
What is NOT kept: the 3D coordinates in the file. dockflow.py always generates its
own 3D start conformer from the SMILES (ETKDG + MMFF), so the input pose is ignored.
Multi-fragment input (salts/mixtures) is rejected visibly, never silently edited.
"""
import os, re

def _libs():
    from openbabel import pybel  # must come before rdkit (swig clash)
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    return pybel, Chem

SUPPORTED = (".smi", ".smiles", ".csv", ".tsv", ".txt", ".mol", ".sdf", ".sd", ".mol2")

def _skeleton(Chem, m):
    """Connectivity-only fingerprint: element + connectivity, no bond orders, charges, H or stereo."""
    r = Chem.RWMol(Chem.RemoveHs(Chem.Mol(m), sanitize=False))
    for bd in r.GetBonds(): bd.SetBondType(Chem.BondType.SINGLE); bd.SetIsAromatic(False)
    for at in r.GetAtoms():
        at.SetFormalCharge(0); at.SetIsAromatic(False); at.SetNoImplicit(True); at.SetNumExplicitHs(0)
        at.SetNumRadicalElectrons(0); at.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    return Chem.MolToSmiles(r, isomericSmiles=False)

def _svg(Chem, mol):
    try:
        from rdkit.Chem.Draw import rdMolDraw2D
        m = Chem.Mol(mol); m.RemoveAllConformers()
        d = rdMolDraw2D.MolDraw2DSVG(220, 150); d.DrawMolecule(m); d.FinishDrawing()
        return d.GetDrawingText()
    except Exception:
        return ""

def _finish(Chem, mol, name, source, rows, fallback_name, kind="file", other=None, other_failed=False):
    """kind 'file' (MOL/SDF/MOL2) or 'smiles' (text). other = second independent reader's mol.
    level: ok | warn (docked only after the user accepts) | reject (never docked)."""
    nm = (name or "").strip() or fallback_name
    def rej(msg):
        rows.append({"name": nm, "smiles": "", "source": source, "ok": False, "level": "reject",
                     "note": msg, "warnings": [], "svg": ""})
    if mol is None: return rej("could not be read or failed sanitisation (invalid valence/bonding)")
    try:
        warns, infos = [], []
        three_d = bool(kind == "file" and mol.GetNumConformers() and mol.GetConformer().Is3D())
        if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
            return rej("unpaired electrons / missing hydrogens after reading (the file looks incomplete or has wrong bond orders)")
        mol = Chem.RemoveHs(mol)
        nfr = len(Chem.GetMolFrags(mol))
        if nfr > 1:
            return rej(f"{nfr} separate fragments (salt, counter-ion or mixture). Not changed automatically: "
                       "remove the extra fragments yourself and re-add the molecule")
        if mol.GetNumHeavyAtoms() < 2: return rej("fewer than 2 heavy atoms")
        if mol.GetNumHeavyAtoms() > 120: return rej("more than 120 heavy atoms (not a small molecule)")
        if other is not None:
            if _skeleton(Chem, mol) != _skeleton(Chem, other):
                return rej("two independent readers (RDKit and Open Babel) disagree on which atoms are bonded; file not trusted")
        elif kind == "file":
            warns.append("second-reader cross-check (Open Babel) could not run on this file; connectivity not cross-checked")
        smi = Chem.MolToSmiles(mol)
        centres = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
        undef = [c for c in centres if c[1] == "?"]
        has_stereo = bool(re.search(r"[@/\\]", smi))
        if "[H]" in smi:
            warns.append("explicit [H] kept because it defines double-bond geometry (e.g. imine/amidine N-H); this is geometry read from the file, not necessarily a real isomer")
        if has_stereo and three_d:
            warns.append("stereo (chiral centres / double-bond geometry) comes only from the 3D coordinates in the file. If the compound is racemic or its stereo is unknown, this docks one arbitrary form")
        if undef:
            (warns if kind == "file" else infos).append(f"{len(undef)} stereocentre(s) undefined; one arbitrary stereoisomer will be docked")
        charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
        if charged:
            net = sum(a.GetFormalCharge() for a in charged)
            (warns if kind == "file" else infos).append(f"charged atoms present (net {net:+d}); taken exactly as written, no pH adjustment")
        level = "warn" if warns else "ok"
        rows.append({"name": nm, "smiles": smi, "source": source, "ok": True, "level": level,
                     "note": "; ".join(infos), "warnings": warns, "svg": _svg(Chem, mol)})
    except Exception as e:
        rej(str(e))

def _smiles_lines(Chem, text, source, rows):
    n = 0
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"): continue
        parts = [p.strip() for p in re.split(r"[,\t;]", line)] if re.search(r"[,\t;]", line) else line.split(None, 1)
        parts = [p for p in parts if p]
        if not parts: continue
        low = [p.lower() for p in parts]
        if n == 0 and low[0] in ("name", "smiles", "compound", "id") and (len(parts) < 2 or low[1] in ("smiles", "name", "structure")):
            continue
        n += 1
        if len(parts) == 1: name, smi = f"{source}_{n}", parts[0]
        else:
            # accept "name,smiles" or "smiles name"
            a, b = parts[0], parts[1]
            ma, mb = Chem.MolFromSmiles(a), Chem.MolFromSmiles(b)
            if ma is not None and mb is None: smi, name = a, b
            elif mb is not None: name, smi = a, b
            else: name, smi = a, b
        _finish(Chem, Chem.MolFromSmiles(smi), name, source, rows, f"{source}_{n}", kind="smiles")
    if n == 0: rows.append({"name": source, "smiles": "", "source": source, "ok": False, "level": "reject", "note": "no molecules found", "warnings": [], "svg": ""})

def _ob_other(pybel, Chem, fmt, text):
    try:
        ob = pybel.readstring(fmt, text)
        m = Chem.MolFromMolBlock(ob.write("mol"), sanitize=False, removeHs=False)
        return m
    except Exception:
        return None

def _mol2_blocks(text):
    parts = re.split(r"(?m)^(?=@<TRIPOS>MOLECULE)", text)
    return [p for p in parts if p.strip().startswith("@<TRIPOS>MOLECULE")]

def parse_text(filename, text):
    """Parse one file (or pasted text if filename is None). Returns list of row dicts."""
    pybel, Chem = _libs()
    rows = []
    base = os.path.splitext(os.path.basename(filename or "pasted"))[0] or "compound"
    ext = os.path.splitext(filename or "")[1].lower()
    src = filename or "pasted text"
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    if ext in (".sdf", ".sd"):
        blocks = [b for b in re.split(r"(?m)^\$\$\$\$[ \t]*$\n?", text) if b.strip()]
        for i, b in enumerate(blocks, 1):
            title = b.split("\n", 1)[0].strip()
            m = Chem.MolFromMolBlock(b, sanitize=True, removeHs=False)
            _finish(Chem, m, title, src, rows, f"{base}_{i}", other=_ob_other(pybel, Chem, "sdf", b + "\n$$$$\n") if m else None)
    elif ext == ".mol":
        title = text.split("\n", 1)[0].strip()
        t0 = text; m = Chem.MolFromMolBlock(t0, sanitize=True, removeHs=False)
        _finish(Chem, m, title, src, rows, base, other=_ob_other(pybel, Chem, "mol", t0) if m else None)
    elif ext == ".mol2":
        blocks = _mol2_blocks(text)
        if not blocks: rows.append({"name": base, "smiles": "", "source": src, "ok": False, "level": "reject", "warnings": [], "svg": "", "note": "no @<TRIPOS>MOLECULE record found"})
        for i, b in enumerate(blocks, 1):
            lines = b.split("\n")
            title = lines[1].strip() if len(lines) > 1 else ""
            m = Chem.MolFromMol2Block(b, sanitize=True, removeHs=False)
            _finish(Chem, m, title, src, rows, f"{base}_{i}", other=_ob_other(pybel, Chem, "mol2", b) if m else None)
    elif ext in (".smi", ".smiles", ".csv", ".tsv", ".txt", ""):
        # pasted text may also be a pasted MOL/SDF/MOL2 block
        if "@<TRIPOS>MOLECULE" in text: return parse_text("pasted.mol2", text)
        if "$$$$" in text or "M  END" in text: return parse_text("pasted.sdf", text)
        _smiles_lines(Chem, text, base if filename else "pasted", rows)
    else:
        rows.append({"name": base, "smiles": "", "source": src, "ok": False, "level": "reject", "warnings": [], "svg": "",
                     "note": f"unsupported file type '{ext}'. Use .smi .csv .txt .mol .sdf .mol2"})
    return rows

def parse_many(items):
    """items: list of (filename_or_None, text). Returns rows with unique names."""
    rows = []
    for fn, tx in items:
        try: rows.extend(parse_text(fn, tx))
        except Exception as e:
            rows.append({"name": fn or "pasted", "smiles": "", "source": fn or "pasted text", "ok": False, "level": "reject", "warnings": [], "svg": "", "note": f"read error: {e}"})
    seen = {}
    for r in rows:
        r["name"] = re.sub(r"[^\w\-. ()+]", "_", r["name"])[:60] or "compound"
        k = r["name"]; seen[k] = seen.get(k, 0) + 1
        if seen[k] > 1: r["name"] = f"{k}_{seen[k]}"
    return rows
