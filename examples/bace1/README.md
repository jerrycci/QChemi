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
以下命令可完整重現論文 Kirsopp et al., Int. J. Quantum Chem., 2022
的「classical CCSD without active space」輸入前處理。

1️⃣ 生成所有 ligand 輸入檔案

執行整批自動流程（包含 ligand.sdf 與 charges.txt）：

python3 batch_generate_ligands_and_charges.py


此步驟會：

從 examples/bace1/data/*.pdb 自動分離 protein 與 ligand；

生成經胍基陽離子化（NH₂–C(=NH⁺)–NH₂）的 ligand.sdf；

以 pdb2pqr 建立蛋白質電荷；

加入 Na⁺ 中和離子；

自動添加一個最終 校正點電荷（確保總電荷 ≈ 0）。

2️⃣ 驗證 ligand 幾何與帶電狀態

以任一配體（例：109-4DK5）進行檢查：

python3 validate_sdf_topology.py ./output/109-4DK5


預期輸出：

Formal charge (SDF) = +1
Has NH2-C(=NH+)-NH2 headgroup? Yes
✅ ligand.sdf checks complete.


這代表：

配體成功帶正電 (+1)；

胍基中心被正確偵測；

幾何未變形（RMSD ≈ 0）。

3️⃣ 驗證 charges.txt 與總電荷中和性

接著確認點電荷與中和狀態：

python3 validate_charges_relaxed.py ./output/109-4DK5


預期輸出：

Note: 17 extra point charges (ions or correction) present.
Total charge from charges.txt = 0.0000 e
✅ charges.txt basic checks complete.


這代表：

17 顆 Na⁺ 及最終校正點電荷已自動加入；

系統總電荷為 0；

所有蛋白電荷成功轉換。

✅ 驗證通過後

每個 ligand 資料夾（如 output/109-4DK5/）將包含：

ligand.sdf — 可直接輸入 DMET / CCSD 模組

protein.pqr — 含溶液相點電荷分布

charges.txt — 以 AMBER10:EHT 對應點電荷
這些輸入即可用於 classical CCSD without active space 方法的能量計算。
