# BACE1 Ligand Dataset (12 Ligands from Hilpert et al. 2013)

This repository provides references and instructions to obtain the **12 BACE1 inhibitor complexes** used in:

> Kirsopp, J. et al. *Quantum Computational Quantification of Protein–Ligand Interactions*,  
> *Int. J. Quantum Chem.* 2022, e26810.  
> (Data sourced from Hilpert et al., *J. Med. Chem.* 2013, 56, 3980–4010.)

---

## 🧬 Dataset Overview

| Ligand ID | PDB ID | Resolution (Å) | Notes |
|:-----------|:--------|:----------------|:------|
| 1b | 3ZLQ | 1.80 | Reference BACE1–inhibitor complex |
| 14d | 4DJW | 1.90 | Hydroxyethylamine inhibitor |
| 66 | 4DJY | 1.95 | Oxazine derivative |
| 67 | 4DJZ | 2.05 | Polar variant of 66 |
| 68 | 4DK0 | 1.95 | Meta-difluorophenyl analog |
| 69 | 4DK1 | 2.00 | Para-trifluoromethyl analog |
| 70 | 4DK2 | 1.90 | Chlorophenyl analog |
| 88 | 4DK3 | 1.70 | High-affinity oxazine ligand |
| 89 | 4DK4 | 1.70 | Ethyl-substituted analog |
| 109 | 4DK5 | 2.05 | Pyridyl-oxazine derivative |
| 110 | 4DK6 | 1.85 | Difluoro-phenyl analog |
| 111 | 4DK7 | 1.80 | Highest pIC₅₀ (~7.9) |

All twelve complexes are publicly available from the **Protein Data Bank (PDB)** and included in the **PDBBind v2020 refined set**.

---

## 📦 Download Instructions

### 1️⃣ From PDB directly
Each entry can be downloaded via:
```
https://www.rcsb.org/structure/<PDB_ID>
```

Example (for 4DK7):
```bash
wget https://files.rcsb.org/download/4DK7.pdb
```

To batch-download all structures:
```bash
for id in 3ZLQ 4DJW 4DJY 4DJZ 4DK0 4DK1 4DK2 4DK3 4DK4 4DK5 4DK6 4DK7; do
  wget "https://files.rcsb.org/download/${id}.pdb"
done
```

### 2️⃣ Extract Ligands (SDF/MOL2)
On each PDB entry page, open **“Ligand Interaction” → Download Chemical Component (SDF/MOL2)**,  
or use RDKit/Open Babel to extract the bound ligand from the downloaded `.pdb`.

### 3️⃣ From PDBBind (refined set)
If you already have the PDBBind v2020 dataset:
```
PDBBind_v2020/refined-set/<PDB_ID>/<PDB_ID>_ligand.sdf
```

---

## 🧩 Example Directory Structure
```
bace1-ligands/
├── data/
│   ├── 3ZLQ.pdb
│   ├── 4DJW.pdb
│   ├── 4DJY.pdb
│   ├── ...
│   └── 4DK7.pdb
├── ligands/
│   ├── 3ZLQ_ligand.sdf
│   ├── ...
│   └── 4DK7_ligand.sdf
└── README.md
```

---

## 🧪 Usage Example (Python + RDKit)
```python
from rdkit import Chem
ligand = Chem.MolFromMolFile("ligands/4DK7_ligand.sdf")
print(ligand.GetNumAtoms())
```

---

## 📖 References
- Hilpert, H. et al. *J. Med. Chem.* 2013, 56, 3980–4010.  
- Kirsopp, J. et al. *Int. J. Quantum Chem.* 2022, e26810.  
- PDBBind v2020 Dataset: <http://www.pdbbind.org.cn/>

---

## 🪪 License
The structures are distributed under the **PDB Usage Policy**.  
This repository only provides public identifiers and retrieval scripts.
