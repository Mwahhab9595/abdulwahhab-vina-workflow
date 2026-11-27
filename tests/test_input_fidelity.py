import os, sys, unittest
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import ligand_input as L
from openbabel import pybel  # noqa: F401  (import order)
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
FX = os.path.join(HERE, "fixtures")

def read(name):
    return L.parse_many([(name, open(os.path.join(FX, name)).read())])

def nostereo(smi):
    m = Chem.MolFromSmiles(smi); return Chem.MolToSmiles(m, isomericSmiles=False)

class T(unittest.TestCase):
    def test_external_identity(self):
        r = read("aspirin_pubchem2244_3d.sdf")[0]
        self.assertEqual(r["level"], "ok")
        self.assertEqual(Chem.CanonSmiles(r["smiles"]), Chem.CanonSmiles("CC(=O)OC1=CC=CC=C1C(=O)O"))
        r = read("ibuprofen_pubchem3672_3d.sdf")[0]
        self.assertEqual(r["level"], "warn")  # stereo invented from 3D of a racemate record
        self.assertEqual(nostereo(r["smiles"]), nostereo("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"))
        r = read("BEN_ideal.sdf")[0]
        self.assertEqual(r["level"], "warn")
        self.assertEqual(rdMolDescriptors.CalcMolFormula(Chem.AddHs(Chem.MolFromSmiles(r["smiles"]))), "C7H8N2")
    def test_levels(self):
        exp = {"rd_sulfonAmide.mol2": "ok", "rd_fusedRing.mol2": "ok", "rd_Sulfonate.mol2": "warn",
               "rd_chargedAmidine.mol2": "warn", "rd_pyridiniumPhenyl.mol2": "warn", "rd_Noxide.mol2": "warn",
               "ob_mol24.mol2": "warn", "ob_culgi_00.mol2": "reject", "synthetic_no_hydrogens.mol2": "reject"}
        for f, lv in exp.items():
            self.assertEqual(read(f)[0]["level"], lv, f)
    def test_rejections(self):
        r = L.parse_many([(None, "NC(=N)c1ccccc1.Cl\nbad,xyz")])
        self.assertEqual([x["level"] for x in r], ["reject", "reject"])
    def test_roundtrip_stereo_preserved(self):
        for smi in ["N[C@@H](C)C(=O)O", "C/C=C/C(=O)O", "C/C=C\\C(=O)O"]:
            m = Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m, randomSeed=3)
            r = L.parse_many([("x.sdf", Chem.MolToMolBlock(m) + "$$$$\n")])[0]
            self.assertEqual(Chem.CanonSmiles(r["smiles"]), Chem.CanonSmiles(smi))
if __name__ == "__main__": unittest.main()
