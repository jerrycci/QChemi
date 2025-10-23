# BACE1 Ligand Preparation Pipeline (Open-Source MOE Replacement)

This project reproduces the **protein–ligand preprocessing workflow** described in  
*Kirsopp et al., Int. J. Quantum Chemistry (2022)*, using **only open-source tools**.  
It automatically generates each ligand’s **`ligand.sdf`** and environment **`charges.txt`**
(AMBER charges + 10 Å water shell equivalent), suitable for subsequent
quantum-chemical fragment and CCSD calculations.

---

## 🧱 1. Environment Setup (Ubuntu + Python venv)

Run the one-click setup script to install all required tools without Conda.

```bash
chmod +x setup_openchem_env_venv.sh
./setup_openchem_env_venv.sh
```

This installs:

| Component | Purpose |
|------------|----------|
| RDKit 2022.x | Ligand SMILES/SDF conversion |
| PyMOL (open-source) | Visualization |
| MDAnalysis + ParmEd | Topology & charge post-processing |
| pdb2pqr             | Protonation  |
| Open Babel 3 | Format conversion (SDF ↔ PDB) |
| BioPython + Pathlib | PDB structure splitting |

Activate the environment:
```bash
source ~/openchem_env/bin/activate
```

---

## ⚙️ 2. Input Folder Structure

```
examples/bace1/
└── data/
    ├── 109-4DK5.pdb
    ├── 110-4DK6.pdb
    ├── …
    └── 89-4DK4.pdb
```

Each `.pdb` should be a full **protein + ligand complex**.

---

## 🚀 3. Run the Batch Generator

```bash
python3 examples/bace1/batch_generate_ligands_and_charges.py
```

Output example:

```
examples/bace1/
├── data/
│   └── 109-4DK5.pdb
├── 109-4DK5/
│   ├── ligand.pdb
│   ├── ligand.sdf
│   ├── protein.pdb
│   ├── protein.pqr
│   └── charges.txt
├── 110-4DK6/
│   ├── ligand.sdf
│   ├── charges.txt
│   └── …
└── …
```

---

## 🧩 4. Script Descriptions

### 🔧 `setup_openchem_env_venv.sh`
Creates `~/openchem_env` venv and installs:
- **pdb2pqr** via APT
- **RDKit 2022.9.5 + NumPy 1.26** (to avoid NumPy 2 ABI issues)
- **Open Babel 3 wheel** for Python
- **PyMOL**, **BioPython**, **MDAnalysis**, **ParmEd**

### 🧬 `batch_generate_ligands_and_charges.py`
- Iterates through `examples/bace1/data/*.pdb`
- Splits each complex into ligand and protein parts
- Generates `ligand.sdf` via RDKit
- Assigns AMBER charges using `pdb2pqr`
- Extracts point charges (x y z q) → `charges.txt`

---

## ✅ 5. Verification & Post-Processing

Check versions:

```bash
python -c "import rdkit; print('RDKit:', rdkit.__version__)"
pdb2pqr --version
obabel -V
```

Preview structures:

```bash
pymol examples/bace1/109-4DK5/ligand.sdf &
```

---

## 📖 Reference
> Kirsopp J., et al. *Quantum computational quantification of protein–ligand interactions*,  
> **Int. J. Quantum Chem.**, 2022. DOI: 10.1002/qua.26976

---

## 🧠 Notes
- Default pH = 7.0 and AMBER force field charges.  
- To add TIP3P water shell or Na⁺ neutralization, extend the script with `tleap`.  
- Works on Ubuntu 20.04 / 22.04 with Python 3.10 or 3.11.

---

## 🧑‍🔬 Example
```bash
source ~/openchem_env/bin/activate
python3 examples/bace1/batch_generate_ligands_and_charges.py
# Output: 12 ligands each with ligand.sdf and charges.txt ready for CCSD calculation
```
