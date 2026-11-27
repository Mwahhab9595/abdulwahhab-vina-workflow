Abdulwahhab Vina Workflow
=========================
A local molecular docking workflow around AutoDock Vina with redocking validation,
multi-criteria pose review, a 3D pose viewer and 2D interaction diagrams.
Your structures and results stay on your computer. The internet is used only to
download PDB files and to check for updates when you press the button.

Important limits (please read)
------------------------------
- Docking scores are estimates, not measured binding energies. Rankings are a guide for human review.
- Validation (redocking) supports but does not prove that results are reliable for your compounds.
- Protein preparation is automatic and approximate: waters, ions and cofactors are removed; protonation
  states and tautomers of ligands are used as written. Check your system by eye.
- The 2D diagram and interaction lists come from simple distance and geometry rules (no angle check for
  hydrogen bonds). They look similar to Discovery Studio diagrams but are not identical to it.
- Wet-lab or further computational confirmation is needed before scientific conclusions.
- This version was built and tested on Linux. Windows packaging has not been tested on a real Windows
  computer by the author of the tool.

How to cite
-----------
Please cite the tools that produced your results:
- AutoDock Vina 1.2: J. Eberhardt, D. Santos-Martins, A. F. Tillack, S. Forli. AutoDock Vina 1.2.0: New Docking
  Methods, Expanded Force Field, and Python Bindings. J. Chem. Inf. Model. 2021, 61, 3891-3898.
- AutoDock Vina: O. Trott, A. J. Olson. J. Comput. Chem. 2010, 31, 455-461.
- Meeko (ligand preparation): Forli Lab, https://github.com/forlilab/Meeko
- RDKit: G. Landrum et al., RDKit: Open-source cheminformatics, https://www.rdkit.org
- Open Babel: N. M. O'Boyle et al. Open Babel: An open chemical toolbox. J. Cheminform. 2011, 3, 33.
- 3Dmol.js (3D viewer): N. Rego, D. Koes. 3Dmol.js: molecular visualization with WebGL. Bioinformatics 2015, 31, 1322-1324.
- Protein Data Bank: H. M. Berman et al. The Protein Data Bank. Nucleic Acids Res. 2000, 28, 235-242.
- NumPy: C. R. Harris et al. Nature 2020, 585, 357-362. SciPy: P. Virtanen et al. Nat. Methods 2020, 17, 261-272.

Licenses of included components
-------------------------------
- This application (source code written for this tool): GPL-2.0, because it uses Open Babel. See LICENSE.
- AutoDock Vina 1.2.7 (vina.exe, run as a separate program): Apache License 2.0 (licenses/AutoDock-Vina-LICENSE.txt).
- Open Babel (openbabel-wheel): GNU GPL version 2. Source: https://github.com/openbabel/openbabel
- Meeko: LGPL-2.1. RDKit: BSD-3-Clause. NumPy, SciPy: BSD-3-Clause.
- pywebview: BSD-3-Clause. pythonnet: MIT. 3Dmol.js: BSD-3-Clause.
- Python runtime and PyInstaller bootloader: PSF License and GPL-2.0 with the PyInstaller exception.
- Microsoft Edge WebView2 runtime (used by the window, not included): provided by Microsoft, preinstalled on current Windows 10/11.
The complete corresponding source of this program is the public repository named in Settings > Updates.
Full license texts are in the licenses folder and in each package's own distribution.
