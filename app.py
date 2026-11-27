#!/usr/bin/env python3
"""Abdulwahhab Vina Workflow entry point.
  app.py                 standalone window (WebView2); falls back to the browser if the window cannot open
  app.py --browser       force the browser
  app.py --no-open       server only
  app.py --engine ...    internal: run the docking engine (used by the app itself)"""
import sys
if len(sys.argv) > 1 and sys.argv[1] == "--engine":
    import dockflow
    dockflow.main(sys.argv[2:])
    sys.exit(0)
if len(sys.argv) > 1 and sys.argv[1] == "--selftest":
    # headless check used by the build: redock 3PTB/BEN with the bundled Vina and print the result
    import os, json, tempfile, csv
    import dockflow, vina_gui
    out = os.path.join(tempfile.mkdtemp(), "st")
    lib = os.path.join(tempfile.mkdtemp(), "lib.csv"); open(lib, "w").write("name,smiles\nBenzamidine,NC(=N)c1ccccc1\n")
    pdb = os.path.join(vina_gui.HERE, "examples", "3PTB.pdb")
    try: dockflow.main(["--pdb", pdb, "--ligand-id", "BEN", "--library", lib, "--out", out, "--exhaustiveness", "8"])
    except SystemExit as e: print("engine exit:", e)
    v = json.load(open(os.path.join(out, "validation.json")))
    print("SELFTEST", "PASS" if v.get("protocol_validated") else "FAIL", "rmsd", v.get("redock_top_rmsd"), "vina", v.get("redock_top_score"))
    sys.exit(0 if v.get("protocol_validated") else 1)
import threading, webbrowser, os, time
import vina_gui

def main():
    srv, port = vina_gui.make_server()
    url = f"http://localhost:{port}"
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    print(f"Abdulwahhab Vina Workflow running at {url}", flush=True)
    if "--no-open" in sys.argv:
        try:
            while True: time.sleep(3600)
        except KeyboardInterrupt: return
    if "--browser" not in sys.argv:
        try:
            import webview
            webview.create_window(f"{vina_gui.version.APP_NAME} {vina_gui.version.VERSION}", url, width=1360, height=860, min_size=(1000, 650))
            webview.start()     # returns when the window is closed
            return
        except Exception as e:
            print("Window could not open, using the browser instead:", e, flush=True)
    webbrowser.open(url)
    try:
        while True: time.sleep(3600)
    except KeyboardInterrupt: pass

if __name__ == "__main__":
    main()
