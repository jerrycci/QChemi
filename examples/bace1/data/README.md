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

# 🧪 MOE 操作手冊：產生 ligand.sdf 與 charges.txt

這份手冊示範如何根據論文 *Kirsopp et al., Int. J. Quantum Chem. (2022)* 的方法，在 **MOE 2020.09** 或更新版中準備兩個必要的檔案：

* **`ligand.sdf`** — 配體量子區域結構
* **`charges.txt`** — 蛋白（含 10 Å 水殼與 Na⁺ 中和）點電荷檔

---

## 一、準備與資料來源

1. **資料來源**：從 PDB 下載 12 個 BACE1 蛋白–配體複合物結構，例如：

   * 4DK5、4DKE、4DKE、4DKE、4DKE... 共 12 個（對應論文 Figure 1）
2. **目的**：

   * 為每個 ligand 建立對應的 `ligand.sdf` 與 `charges.txt`
   * 確保每個 ligand 都有對應的局部蛋白電荷環境與配體結構

---

## 二、在 MOE 產生 ligand.sdf

### 1️⃣ 載入蛋白–配體複合物

* 開啟 **MOE** → `File > Open` → 選取對應的 PDB (例如 4DK5.pdb)
* 檢查配體是否標示為 `HETATM`

### 2️⃣ 準備配體結構

* 選取配體：`Select > Ligand > Current Site`
* 執行：`Structure > 3D Protonate`

  * 使用預設設定 (pH 7.0)
  * 確保所有氫原子已加上

### 3️⃣ 確認配體電荷

* 執行：`Compute > Partial Charges`

  * 選擇 Force Field: **AMBER10:EHT**
  * 檢查配體總電荷是否為 +1（因 NH₂–C(=NH⁺)–NH₂ 結構）

### 4️⃣ 匯出 ligand.sdf

* 選擇：`File > Export > Molecule Format = MDL SD File (*.sdf)`
* 命名為 `ligands12_real/<ligid>.sdf`，例如：`ligands12_real/109.sdf`

---

## 三、在 MOE 產生 charges.txt

### 1️⃣ 加入氫原子與水殼

* `Compute > 3D Protonate` → 預設設定。
* `Applications > Structure Preparation > Solvate`

  * **Shell Radius = 10 Å**
  * **Solvent = TIP3P Water**
* MOE 會自動生成水分子與 Na⁺ 離子中和系統電荷。

### 2️⃣ 指派 AMBER10:EHT 電荷

* 執行：`Compute > Partial Charges`

  * **Force Field:** AMBER10:EHT
  * **Assign Charges To:** All Atoms (protein + water + Na⁺)

### 3️⃣ 匯出電荷表

* 開啟：`Database Viewer`
* `File > Export > CSV`
* 選取欄位：**Atom Index, x, y, z, Charge**
* 儲存為：`charges.csv`

### 4️⃣ 轉換成 charges.txt

在 Python 執行以下腳本：

```python
import pandas as pd

df = pd.read_csv("charges.csv")
df[['x', 'y', 'z', 'Charge']].to_csv("charges.txt", sep=' ', index=False, header=False)
print(f"Wrote {len(df)} charges to charges.txt")
```

產出格式：

```
12.435  15.218  8.433  0.123
12.912  16.041  7.812 -0.234
...
```

此檔案即為 PySCF `qmmm.mm_charge()` 所需的輸入。

---

## 四、檔案結構建議

```
examples/bace1/
├── data/
│   ├── 109-4DK5.pdb
│   ├── 110-4DKE.pdb
│   └── ...
├── ligands12_real/
│   ├── 109.sdf
│   ├── 110.sdf
│   └── ...
├── 109-4DK5/
│   ├── protein/
│   │   └── charges.txt
│   └── solvent/
└── 110-4DKE/
    ├── protein/
    │   └── charges.txt
    └── solvent/
```

---

## 五、批次轉換工具 convert_moe_charges.py

這個工具可批次將多個 ligand 的 MOE 匯出 CSV 轉成對應的 charges.txt。

```python
#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import argparse

def convert_one(csv_path: Path):
    df = pd.read_csv(csv_path)
    txt_path = csv_path.with_suffix('.txt')
    df[['x','y','z','Charge']].to_csv(txt_path, sep=' ', index=False, header=False)
    print(f"[OK] {csv_path.name} → {txt_path.name} ({len(df)} atoms)")

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--indir', required=True, help='Folder containing MOE-exported charges.csv files')
    args = p.parse_args()
    indir = Path(args.indir)

    for csv in indir.glob('*.csv'):
        convert_one(csv)

if __name__ == '__main__':
    main()
```

使用方式：

```bash
python3 convert_moe_charges.py --indir ./examples/bace1/charges_raw
```

轉換後會自動產生 `*.txt` 給每個 ligand 使用。

---

## 六、確認與建議

* 每個 ligand 都必須擁有：

  * `ligand.sdf`（3D Protonate 後導出）
  * `charges.txt`（AMBER10:EHT + TIP3P 水殼 + Na⁺ 中和）
* 這兩者組成了論文中所需的：

  * ( E_{ligand-in-protein(aq)} ) 與 ( E_{ligand-in-solvent(aq)} ) 的輸入。

