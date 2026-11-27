# What the flow does (read from the code, dockflow.py unchanged)

1. Load the target: PDB ID (downloaded from RCSB) or an uploaded .pdb.
2. Chains (GUI step, before dockflow): you choose the chain(s); other chains are removed from the file handed to dockflow. dockflow itself keeps every chain it receives.
2b. Receptor preparation: keep ATOM records only (alternate location blank/A). All HETATM (waters, metals, cofactors, ions, bound ligand) are removed. All chains are kept. Open Babel adds hydrogens (pH 7.4 correction) and writes a rigid receptor PDBQT. NOT done: repair of missing residues/atoms, per-residue protonation or histidine tautomer checks, minimisation, keeping essential waters/metals/cofactors, flexible residues. Review these by hand.
3. Pocket and box: the bound ligand (largest non-solvent HET, or --ligand-id) + 5 A padding (20-40 A box; errors instead of clipping). Without a bound ligand: reference SMILES plus a geometric cavity finder (verify manually).
4. Validation gate: rebuild the bound ligand from the RCSB ideal SDF, re-dock it from a random conformer, pass if top-pose RMSD <= 2.0 A. Otherwise stop (unless exploratory mode, which labels outputs UNVALIDATED).
5. Ligand preparation, per compound, from SMILES: RDKit add hydrogens, 3D embedding (ETKDGv3), MMFF optimisation (one start conformer), Meeko (atom types, Gasteiger charges, rotatable bonds) to PDBQT. NOT done: pH protonation/ionisation (docked as written), tautomer or stereoisomer enumeration, salt removal, multiple conformers.
6. Docking with Vina scoring (exhaustiveness 8, 10 poses by default).
7. Pose review: Vina score, ligand efficiency, key-residue contacts, N/O polar contacts (distance only), pose consensus, MMFF strain, burial, clashes; weighted composite (unbenchmarked defaults). Outputs: ranked.csv/json, all_poses_audit.csv, best_poses.pdb, complexes, validation.json, report.md, plus input_check.csv from the GUI.
The ranking is a heuristic shortlist, not proof of activity.
