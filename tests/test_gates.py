"""Gate tests (need internet for RCSB; ~2 min total). Run: python -m unittest tests/test_gates.py -v"""
import os, subprocess, sys, tempfile, unittest, json
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(HERE, "dockflow.py")
LIB = os.path.join(HERE, "trypsin_apo_lib.csv")

def run(*args):
    out = tempfile.mkdtemp()
    p = subprocess.run([sys.executable, SCRIPT, "--library", LIB, "--out", out, *args],
                       capture_output=True, text=True)
    return p, out

class Gates(unittest.TestCase):
    def test_validation_pass_produces_outputs(self):
        p, out = run("--pdb", "3PTB", "--ligand-id", "BEN")
        self.assertEqual(p.returncode, 0, p.stderr[-500:])
        v = json.load(open(f"{out}/validation.json"))
        self.assertTrue(v["protocol_validated"]); self.assertLessEqual(v["redock_top_rmsd"], 2.0)
        for f in ("ranked.csv", "best_poses.pdb", "all_poses_audit.csv"):
            self.assertTrue(os.path.exists(f"{out}/{f}"), f)
        self.assertFalse(os.path.exists(f"{out}/UNVALIDATED.txt"))

    def test_failed_validation_blocks_screening(self):
        # impossible cutoff forces a validation failure
        p, out = run("--pdb", "3PTB", "--ligand-id", "BEN", "--rmsd-cutoff", "0.01")
        self.assertNotEqual(p.returncode, 0)
        self.assertFalse(os.path.exists(f"{out}/ranked.csv"))
        self.assertFalse(os.path.exists(f"{out}/best_poses.pdb"))
        self.assertEqual(json.load(open(f"{out}/validation.json"))["status"], "UNVALIDATED")

    def test_override_is_labelled_unvalidated(self):
        p, out = run("--pdb", "3PTB", "--ligand-id", "BEN", "--rmsd-cutoff", "0.01", "--allow-unvalidated")
        self.assertEqual(p.returncode, 0, p.stderr[-500:])
        self.assertTrue(os.path.exists(f"{out}/UNVALIDATED.txt"))
        self.assertIn("UNVALIDATED", open(f"{out}/report.md").read().splitlines()[0])
        self.assertIn("UNVALIDATED", open(f"{out}/best_poses.pdb").readline())

    def test_apo_without_override_stops(self):
        # no co-crystal pose: cannot be validated, must stop by default
        pdb = os.path.join(tempfile.mkdtemp(), "apo.pdb")
        import urllib.request
        txt = urllib.request.urlopen("https://files.rcsb.org/download/3PTB.pdb").read().decode()
        open(pdb, "w").write("\n".join(l for l in txt.splitlines() if not l.startswith("HETATM")))
        p, out = run("--pdb", pdb, "--ref-smiles", "NC(=N)c1ccccc1")
        self.assertNotEqual(p.returncode, 0)
        self.assertFalse(os.path.exists(f"{out}/ranked.csv"))

    def test_box_never_silently_clipped(self):
        p, out = run("--pdb", "3PTB", "--ligand-id", "BEN", "--max-box", "10")
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("exceeds --max-box", p.stdout + p.stderr)
        self.assertFalse(os.path.exists(f"{out}/ranked.csv"))

    def test_reused_output_dir_rejected(self):
        out = tempfile.mkdtemp()
        open(os.path.join(out, "ranked.csv"), "w").write("stale")
        p = subprocess.run([sys.executable, SCRIPT, "--pdb", "3PTB", "--ligand-id", "BEN", "--library", LIB,
                            "--out", out, "--rmsd-cutoff", "0.01"], capture_output=True, text=True)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("not empty", p.stdout + p.stderr)
        self.assertEqual(open(os.path.join(out, "ranked.csv")).read(), "stale")  # untouched, not trusted

if __name__ == "__main__":
    unittest.main()
