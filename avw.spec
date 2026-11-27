# PyInstaller spec. Usage: pyinstaller avw.spec            (folder build)
#                         AVW_ONEFILE=1 pyinstaller avw.spec   (single .exe)
import os
from PyInstaller.utils.hooks import collect_all
datas = [("gui", "gui"), ("examples", "examples"), ("bin", "bin"), ("licenses", "licenses"),
         ("ABOUT.md", "."), ("LICENSE", "."), ("trypsin_library.csv", ".")]
binaries = []; hidden = ["dockflow", "avw_analysis", "ligand_input", "vina_backend", "vina_gui", "version", "webview"]
for pkg in ("openbabel", "meeko", "rdkit", "scipy", "webview"):
    try:
        d, b, h = collect_all(pkg); datas += d; binaries += b; hidden += h
    except Exception as e:
        print("collect_all failed for", pkg, e)
a = Analysis(["app.py"], pathex=["."], binaries=binaries, datas=datas, hiddenimports=hidden, excludes=["tkinter", "matplotlib"])
pyz = PYZ(a.pure)
name = "AbdulwahhabVinaWorkflow"
if os.environ.get("AVW_ONEFILE"):
    exe = EXE(pyz, a.scripts, a.binaries, a.datas, name=name, console=False, upx=False)
else:
    exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name=name, console=False, upx=False)
    coll = COLLECT(exe, a.binaries, a.datas, name=name, upx=False)
