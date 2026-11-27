"""Structure analysis helpers for Abdulwahhab Vina Workflow:
chain suggestions, pocket suggestions, box preview, and 2D interaction diagrams."""
import math, re
import numpy as np
from scipy.spatial import cKDTree

AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
       "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
HYDROPHOBIC = {"ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO","CYS","TYR"}
POS = {"LYS","ARG","HIS"}; NEG = {"ASP","GLU"}
AROM = {"PHE":["CG","CD1","CD2","CE1","CE2","CZ"], "TYR":["CG","CD1","CD2","CE1","CE2","CZ"],
        "TRP":["CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"], "HIS":["CG","ND1","CD2","CE1","NE2"]}
NOT_LIGAND = {"HOH","WAT","DOD","SO4","PO4","GOL","EDO","PEG","PGE","MPD","DMS","ACT","ACE","NAG","BMA","MAN","FMT","TRS",
              "EPE","MES","CL","NA","K","MG","ZN","CA","MN","FE","CU","NI","CO","CD","IOD","BR","NO3","SCN","IMD","BME",
              "CIT","TAR","MLI","PG4","1PE"}

def _xyz(l): return [float(l[30:38]), float(l[38:46]), float(l[46:54])]

def chain_sequences(text):
    seq = {}; seen = set()
    for l in text.splitlines():
        if l.startswith("ENDMDL"): break
        if l.startswith("ATOM") and l[16] in " A" and l[12:16].strip() == "CA":
            key = (l[21], l[22:27])
            if key in seen: continue
            seen.add(key); seq.setdefault(l[21], []).append(AA3.get(l[17:20], "X"))
    return {c: "".join(v) for c, v in seq.items()}

def suggest_chains(info, text):
    """Rule-based chain suggestion. info = inspect result (chains, ligands). Returns dict(chains, reason, groups)."""
    seqs = chain_sequences(text); ids = [c["id"] for c in info["chains"]]
    groups = {}
    for c in ids: groups.setdefault(seqs.get(c, c), []).append(c)
    groups = list(groups.values())
    if len(ids) == 1:
        return {"chains": ids, "reason": f"Only one protein chain ({ids[0]}) in this structure.", "groups": groups}
    ligs = info.get("ligands") or []
    if ligs:
        lg = ligs[0]; tot = sum(lg["contacts"].values()) or 1
        share = {c: n / tot for c, n in lg["contacts"].items()}
        main = [c for c, s in share.items() if s >= 0.30]
        if len(main) >= 2:
            return {"chains": sorted(main), "groups": groups, "reason":
                    f"Ligand {lg['res']} contacts chains {', '.join(sorted(main))} with at least 30% of its contact residues each, so its site is shared between them. Both are suggested. Check this in the 3D view."}
        if main:
            c = main[0]
            same = [g for g in groups if c in g and len(g) > 1]
            extra = f" Chains {', '.join(same[0])} have identical sequences (copies); the one with the ligand was chosen." if same else ""
            return {"chains": [c], "groups": groups, "reason":
                    f"Ligand {lg['res']} sits mostly in chain {c} ({round(share[c]*100)}% of its contact residues within 4 A).{extra}"}
    for g in groups:
        if len(g) > 1 and len(groups) == 1:
            return {"chains": [g[0]], "groups": groups, "reason":
                    f"Chains {', '.join(g)} have identical sequences, so they are copies. The first ({g[0]}) is suggested. If the active site is built from two chains (dimer), select both."}
    big = max(info["chains"], key=lambda c: c["residues"])
    return {"chains": [big["id"]], "groups": groups, "reason":
            f"No ligand found and the chains differ in sequence. The largest chain ({big['id']}, {big['residues']} residues) is suggested. Please check which chain holds your site of interest."}

def pocket_suggestions(prot_lines, find_pockets, lining_residues, top=3):
    pk = find_pockets(np.array([_xyz(l) for l in prot_lines]), top=top)
    out = []
    for i, p in enumerate(pk):
        pts = p["points"]; ext = pts.max(0) - pts.min(0)
        res = lining_residues(prot_lines, pts, 4.5)
        out.append({"rank": i + 1, "points": int(p["n_points"]), "volume_A3": int(round(p["n_points"] * 1.5 ** 3)),
                    "center": [round(float(x), 2) for x in p["center"]],
                    "extent": [round(float(x), 1) for x in ext],
                    "residues": res[:40], "n_residues": len(res)})
    return out

# ------------------------------------------------------------------ interactions
def _ring_info(coords):
    c = coords.mean(0); u, s, vt = np.linalg.svd(coords - c); return c, vt[2]

def detect_interactions(mol, prot_lines):
    """mol: RDKit heavy-atom mol with a conformer; prot_lines: PDB ATOM lines. Returns list of dicts."""
    pos = mol.GetConformer().GetPositions()
    P = np.array([_xyz(l) for l in prot_lines]); tree = cKDTree(P)
    res_of = [(l[21], l[22:26].strip(), l[17:20]) for l in prot_lines]
    names = [l[12:16].strip() for l in prot_lines]; el = [(l[76:78].strip() or l[12:16].strip()[0]).upper() for l in prot_lines]
    out = []; seen = set()
    def add(kind, ai, rk, d, pp=None, lp=None):
        k = (kind, ai, rk)
        if k in seen: return
        seen.add(k); out.append({"type": kind, "lig_atom": int(ai), "res": f"{rk[2]}{rk[1]}", "chain": rk[0],
                                  "resname": rk[2], "dist": round(float(d), 2),
                                  "lpos": [round(float(x), 3) for x in (pos[ai] if lp is None else lp)],
                                  "ppos": [round(float(x), 3) for x in pp]})
    charged_lig = {}
    for a in mol.GetAtoms():
        i = a.GetIdx(); s = a.GetSymbol()
        if s in ("N", "O"):
            for j in tree.query_ball_point(pos[i], 3.5):
                if el[j] in ("N", "O"):
                    d = np.linalg.norm(P[j] - pos[i])
                    if d >= 2.4: add("hbond", i, res_of[j], d, P[j])
        if s == "C":
            for j in tree.query_ball_point(pos[i], 4.0):
                if el[j] == "C" and res_of[j][2] in HYDROPHOBIC and names[j] not in ("C", "CA", "O", "N"):
                    d = np.linalg.norm(P[j] - pos[i]); add("hydrophobic", i, res_of[j], d, P[j])
        # ionic: ligand N+ or acidic O with opposite protein charge
        pos_lig = (s == "N" and (a.GetFormalCharge() > 0 or any(n.GetSymbol() == "C" and any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(n).GetSymbol() == "N" for b in n.GetBonds()) for n in a.GetNeighbors()) and not a.GetIsAromatic()))
        neg_lig = (s == "O" and a.GetFormalCharge() < 0) or (s == "O" and any(n.GetSymbol() in ("C", "S", "P") and sum(1 for m in n.GetNeighbors() if m.GetSymbol() == "O") >= 2 for n in a.GetNeighbors()) and a.GetTotalNumHs() == 0 and a.GetDegree() == 1)
        for j in tree.query_ball_point(pos[i], 4.0):
            rn = res_of[j][2]; d = np.linalg.norm(P[j] - pos[i])
            if pos_lig and rn in NEG and el[j] == "O" and names[j] not in ("O", "OXT"): add("ionic", i, res_of[j], d, P[j])
            if neg_lig and rn in POS and el[j] == "N" and names[j] not in ("N",): add("ionic", i, res_of[j], d, P[j])
    # pi stacking
    ri = mol.GetRingInfo()
    prot_rings = {}
    for idx, l in enumerate(prot_lines):
        rn = l[17:20]
        if rn in AROM and names[idx] in AROM[rn]: prot_rings.setdefault((l[21], l[22:26].strip(), rn), []).append(P[idx])
    for ring in ri.AtomRings():
        if not all(mol.GetAtomWithIdx(k).GetIsAromatic() for k in ring): continue
        c1, n1 = _ring_info(pos[list(ring)])
        for rk, pts in prot_rings.items():
            if len(pts) < 5: continue
            c2, n2 = _ring_info(np.array(pts)); d = np.linalg.norm(c1 - c2)
            if d > 5.5: continue
            ang = math.degrees(math.acos(min(1, abs(float(np.dot(n1, n2))))))
            if ang < 30 or ang > 60: add("pi", ring[0], rk, d, c2, c1)
    return out

COLORS = {"hbond": "#2e9d4f", "hydrophobic": "#e38b1d", "pi": "#8a4fd6", "ionic": "#d8453d"}
LABELS = {"hbond": "Hydrogen bond", "hydrophobic": "Hydrophobic", "pi": "Pi stacking", "ionic": "Ionic / salt bridge"}
def _rescolor(rn):
    if rn in HYDROPHOBIC: return "#bfe6b0"
    if rn in POS: return "#a9c4f5"
    if rn in NEG: return "#f5b3b0"
    return "#d9c2ee"

def diagram_svg(mol, prot_lines, title=""):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
    inter = detect_interactions(mol, prot_lines)
    m2 = Chem.Mol(mol); m2.RemoveAllConformers(); AllChem.Compute2DCoords(m2)
    W, H = 900, 760
    d = rdMolDraw2D.MolDraw2DSVG(W, H)
    d.drawOptions().padding = 0.32; d.drawOptions().addAtomIndices = False
    d.DrawMolecule(m2); d.FinishDrawing()
    svg = d.GetDrawingText()
    px = {i: d.GetDrawCoords(i) for i in range(m2.GetNumAtoms())}
    cen = np.array([np.mean([px[i].x for i in px]), np.mean([px[i].y for i in px])])
    # group interactions by residue
    byres = {}
    for it in inter: byres.setdefault((it["chain"], it["res"]), []).append(it)
    items = []
    for (ch, rs), its in byres.items():
        pts = np.array([[px[i["lig_atom"]].x, px[i["lig_atom"]].y] for i in its]).mean(0)
        v = pts - cen; ang = math.atan2(v[1], v[0]) if np.linalg.norm(v) > 1e-6 else 0.0
        items.append([ang, ch, rs, its])
    items.sort(key=lambda t: t[0])
    n = len(items); R = np.array([W * 0.40, H * 0.40])
    angs = [t[0] for t in items]
    minsep = min(0.62, 2 * math.pi / max(n, 1))
    if n > 1:
        mean0 = sum(angs) / n
        for k in range(1, n):
            if angs[k] - angs[k - 1] < minsep: angs[k] = angs[k - 1] + minsep
        if angs[-1] - angs[0] > 2 * math.pi - minsep:       # too crowded: spread evenly
            angs = [angs[0] + k * 2 * math.pi / n for k in range(n)]
        shift = mean0 - sum(angs) / n
        angs = [a + shift for a in angs]
    extra = []
    for k, (a0, ch, rs, its) in enumerate(items):
        a = angs[k]; cx = W / 2 + R[0] * math.cos(a); cy = H / 2 + R[1] * math.sin(a)
        kinds = sorted({i["type"] for i in its}, key=lambda t: list(COLORS).index(t))
        rn = its[0]["resname"]; label = rn.capitalize() + rs[len(rn):] + (f":{ch}" if ch.strip() else "")
        for ki, kind in enumerate(kinds):
            ii = [i for i in its if i["type"] == kind]; best = min(ii, key=lambda i: i["dist"])
            p = px[best["lig_atom"]]
            vx, vy = cx - p.x, cy - p.y; ln = math.hypot(vx, vy) or 1.0
            ox, oy = -vy / ln * 3.5 * (ki - (len(kinds) - 1) / 2), vx / ln * 3.5 * (ki - (len(kinds) - 1) / 2)
            dash = ' stroke-dasharray="6,4"' if kind in ("hbond", "pi") else ' stroke-dasharray="2,3"' if kind == "hydrophobic" else ""
            extra.append(f'<line x1="{p.x+ox:.1f}" y1="{p.y+oy:.1f}" x2="{cx+ox:.1f}" y2="{cy+oy:.1f}" stroke="{COLORS[kind]}" stroke-width="2.2"{dash} opacity="0.9"/>')
            if kind == "hbond":
                mx, my = (p.x + cx) / 2, (p.y + cy) / 2
                extra.append(f'<text x="{mx:.1f}" y="{my-3:.1f}" font-size="11" fill="{COLORS[kind]}" text-anchor="middle" font-family="Arial">{best["dist"]:.1f}</text>')
        extra.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="31" fill="{_rescolor(rn)}" stroke="#555" stroke-width="1.2"/>')
        extra.append(f'<text x="{cx:.1f}" y="{cy+4:.1f}" font-size="12.5" font-weight="bold" text-anchor="middle" fill="#222" font-family="Arial">{label}</text>')
    # legend
    leg = []
    y0 = 24
    leg.append(f'<text x="14" y="{y0}" font-size="14" font-weight="bold" font-family="Arial" fill="#222">{title}</text>')
    present = [k for k in COLORS if any(i["type"] == k for i in inter)]
    for q, k in enumerate(present):
        yy = H - 18 - 20 * (len(present) - 1 - q)
        leg.append(f'<line x1="14" y1="{yy-4}" x2="48" y2="{yy-4}" stroke="{COLORS[k]}" stroke-width="2.4"/><text x="56" y="{yy}" font-size="12" font-family="Arial" fill="#222">{LABELS[k]}</text>')
    body = "".join(extra + leg)
    svg = svg.replace("</svg>", body + "</svg>")
    return svg, inter
