# What the input reader keeps and what it does not (tested)

Pipeline: file -> RDKit parse (sanitised) -> remove H -> canonical SMILES -> dockflow builds its own 3D from the SMILES. Tests: tests/test_input_fidelity.py.

| Property | Preserved? | Evidence / limit |
|---|---|---|
| Atoms and connectivity | Yes, cross-checked by a second reader (Open Babel) on MOL/SDF/MOL2; mismatch = reject | Bond orders and charges are NOT cross-checked: Open Babel's re-export lost charges/aromaticity on charged MOL2 files, so only connectivity is compared |
| Bond orders, aromaticity | As written; RDKit must sanitise or the entry is rejected | culgi_00.mol2 (real file) is rejected |
| Formal charges | As written, no pH change | any charged atom in a file = warning you must accept. A charged Corina MOL2 is read as charged |
| Hydrogens | Removed, then re-added by dockflow (explicit H kept only when it defines double-bond geometry, flagged) | A MOL2 with no hydrogens gives unpaired electrons and is rejected |
| Chirality | Round trip of a 3D file made from a SMILES with defined stereo returns the same stereo (tested: L-alanine, E and Z butenoate). Flagged anyway for 3D-only stereo | PubChem's 3D ibuprofen record (a racemate) reads as one enantiomer: the reader cannot know it was a racemate. Flagged as warning |
| Double-bond E/Z | Same as chirality; imine/amidine N-H geometry from coordinates is flagged | RCSB BEN (benzamidine) reads as [H]/N=C(/N)... : flagged |
| Undefined stereocentres | Not enumerated; one arbitrary stereoisomer is docked | flagged for files, noted for SMILES text |
| Conformation (3D coordinates) of the file | NOT used | dockflow rebuilds a conformer from the SMILES |
| Charge state, protonation and tautomer | USED exactly as written in the file | nothing is enumerated or pH-adjusted; a different tautomer/ionisation state must be supplied by you |
| Salts, mixtures | Rejected, never edited | |
| Names | Taken from the SDF title / MOL2 name / CSV column; duplicates get _2 suffixes | |
| Metals, >120 heavy atoms | Not supported (>120 heavy atoms rejected; metals untested) | |
