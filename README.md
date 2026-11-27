# Abdulwahhab Vina Workflow

Local AutoDock Vina docking workflow: protein and ligand preparation, redocking validation (stops if it fails),
multi-criteria pose review, 3D pose viewer, 2D interaction diagram, Excel/CSV, PDB and image export.

## Windows (students)
Download the release zip, extract it, double-click `AbdulwahhabVinaWorkflow.exe`. Nothing to install
(Vina is bundled). A standalone window opens; if it cannot, your browser opens instead.
Status: built by GitHub Actions; see the release notes for whether it has been tested on real Windows.

## From source (Linux / macOS / Windows)
    pip install -r requirements-build.txt
    python app.py               # window; add --browser to use the browser
    python app.py --selftest    # redock 3PTB/BEN and print PASS/FAIL
Set `AVW_VINA` to a Vina executable, or put it in `bin/` (`bin/vina_1.2.7_win.exe` on Windows, `bin/vina` elsewhere).

## Build
`.github/workflows/build-windows.yml` (manual run only) builds the folder version and a single-file exe on a Windows runner and
runs `--selftest` on both before publishing.

## Updates
Settings > Check for updates reads the latest release of the repository named in `version.py`. Nothing installs automatically.

## License
GPL-2.0 (see LICENSE) because Open Babel (GPL-2.0) is used. Third-party licenses: `ABOUT.md` and `licenses/`.
Docking scores are estimates; results need human review and experimental confirmation.
