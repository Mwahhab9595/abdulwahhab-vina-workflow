"""Vina backend: runs the AutoDock Vina command-line program (vina.exe on Windows).

The program is searched in this order:
  1. environment variable AVW_VINA (a full path)
  2. a bundled copy: bin/vina_1.2.7_win.exe (Windows) or bin/vina (Linux/macOS)
     next to the app, or inside the PyInstaller bundle
  3. 'vina' on PATH
"""
import os, shutil, subprocess, sys, tempfile, re

def _roots():
    r = []
    if getattr(sys, "_MEIPASS", None): r.append(sys._MEIPASS)
    r.append(os.path.dirname(os.path.abspath(__file__)))
    r.append(os.path.dirname(sys.executable))
    return r

def find_vina():
    p = os.environ.get("AVW_VINA")
    if p and os.path.isfile(p): return p
    names = ["bin/vina_1.2.7_win.exe", "bin/vina.exe"] if os.name == "nt" else ["bin/vina", "bin/vina_linux"]
    for root in _roots():
        for n in names:
            f = os.path.join(root, n)
            if os.path.isfile(f): return f
    return shutil.which("vina")

def vina_version(exe=None):
    exe = exe or find_vina()
    if not exe: return None
    try:
        r = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=20, **_flags())
        return (r.stdout or r.stderr).strip().splitlines()[0]
    except Exception:
        return None

def _flags():
    # no console window flashing on Windows
    return {"creationflags": 0x08000000} if os.name == "nt" else {}

def run_vina_cli(rec_qt, lig_qt_text, center, size, exh=8, npose=10, seed=42, cpu=0,
                 scoring="vina", energy_range=3.0, min_rmsd=1.0, exe=None, spacing=0.375):
    """Dock one ligand. Returns (output_pdbqt_text, [energies])."""
    exe = exe or find_vina()
    if not exe:
        raise RuntimeError("Vina program not found. Reinstall the app or set the Vina path in Settings.")
    with tempfile.TemporaryDirectory(prefix="avw_") as td:
        lq = os.path.join(td, "lig.pdbqt"); oq = os.path.join(td, "out.pdbqt")
        open(lq, "w").write(lig_qt_text)
        cmd = [exe, "--receptor", os.path.abspath(rec_qt), "--ligand", lq,
               "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}", "--center_z", f"{center[2]:.3f}",
               "--size_x", f"{size[0]:.3f}", "--size_y", f"{size[1]:.3f}", "--size_z", f"{size[2]:.3f}",
               "--exhaustiveness", str(int(exh)), "--num_modes", str(int(npose)), "--seed", str(int(seed)),
               "--cpu", str(int(cpu)), "--scoring", scoring, "--energy_range", str(energy_range),
               "--min_rmsd", str(min_rmsd), "--spacing", str(spacing), "--out", oq]
        r = subprocess.run(cmd, capture_output=True, text=True, **_flags())
        if r.returncode != 0 or not os.path.exists(oq):
            raise RuntimeError("Vina failed: " + (r.stderr or r.stdout).strip()[-400:])
        txt = open(oq).read()
    en = [float(m.group(1)) for m in re.finditer(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", txt)]
    if not en: raise RuntimeError("Vina produced no poses (try a bigger box)")
    return txt, en
