# Input-fidelity fixtures: where each file comes from and what it should read as

Downloaded 2026-11-25. "External" = file written by someone else's software, not by this project.

| File | Source (URL) | Expected identity and how it is compared | Expected GUI level |
|---|---|---|---|
| BEN_ideal.sdf | RCSB ideal coordinates for ligand BEN: https://files.rcsb.org/ligands/download/BEN_ideal.sdf | Benzamidine, formula C7H8N2 (RCSB chemical component record https://data.rcsb.org/rest/v1/core/chemcomp/BEN). Compared by formula and heavy-atom skeleton | warn (explicit [H]/C=N geometry taken from coordinates) |
| aspirin_pubchem2244_3d.sdf | PubChem CID 2244, 3D record: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/SDF?record_type=3d | PubChem SMILES CC(=O)OC1=CC=CC=C1C(=O)O (https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/IsomericSMILES/JSON). Compared by canonical SMILES | ok |
| ibuprofen_pubchem3672_3d.sdf | PubChem CID 3672, 3D record | PubChem SMILES CC(C)CC1=CC=C(C=C1)C(C)C(=O)O has NO stereo (racemate record). The 3D file still contains coordinates, so the reader invents one enantiomer. Compared ignoring stereo; stereo flagged | warn |
| rd_sulfonAmide.mol2, rd_fusedRing.mol2, rd_Sulfonate.mol2, rd_chargedAmidine.mol2, rd_pyridiniumPhenyl.mol2, rd_Noxide.mol2 | RDKit test data (Corina-written MOL2): https://github.com/rdkit/rdkit/tree/master/Code/GraphMol/FileParsers/test_data | Expected identity is OUR reading of the file's atoms (methanesulfonamide, benzo[g]quinoxaline, methanesulfonate anion, 1,3-dimethylbenzimidazolium, N-phenylpyridinium, an azoxy compound), cross-checked only by a second reader (Open Babel) for connectivity. Not externally verified | ok / warn as in test |
| ob_mol24.mol2 | Open Babel test files: https://github.com/openbabel/openbabel/tree/master/test/files | Has one stereocentre defined only by 3D coordinates | warn |
| ob_culgi_00.mol2 | Open Babel test files | RDKit cannot sanitise it (valence error) | reject |
| synthetic_no_hydrogens.mol2 | made by this project | deliberately incomplete MOL2 without hydrogens | reject |
