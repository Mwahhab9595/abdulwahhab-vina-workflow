# Sources (pages actually opened/returned by search, Nov 25 2026)

Software
- AutoDock Vina 1.2.0 (Eberhardt et al., J Chem Inf Model 2021): https://pubs.acs.org/doi/10.1021/acs.jcim.1c00203 ; PubMed https://pubmed.ncbi.nlm.nih.gov/34278794
- AutoDock Vina code/docs: https://github.com/ccsb-scripps/AutoDock-Vina
- Meeko (ligand preparation to PDBQT): https://github.com/forlilab/Meeko
- RDKit and Open Babel are used via the pip packages `rdkit` and `openbabel-wheel` (no page opened for them).

Protocol
- Redocking the co-crystallised ligand to validate a protocol, as a step that is often neglected: https://pmc.ncbi.nlm.nih.gov/articles/PMC12064015/ (Ten quick tips for meaningful, reproducible docking). Text fetched; it states the principle, not a numeric cutoff.
- The 2 A RMSD cutoff is a widely used convention. Only a search snippet was seen for this, not verified in the page text: https://proteiniq.io/guides/redocking-docking-validation

Experimental test case
- 3PTB, trypsin with benzamidine (co-crystal ligand BEN, used for the redocking test): https://www.rcsb.org/structure/3PTB ; paper: https://journals.iucr.org/paper?a22146=

Not implemented: automated literature lookup of reference inhibitors/activators. The user supplies reference SMILES.
