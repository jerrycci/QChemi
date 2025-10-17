
# protein-ligand-qiskit-dmet (V2 primitives, DMET scaffold)

本範例以 **Qiskit 2.x**/**Qiskit Nature 0.7** 的 **V2 primitives** 路線，提供：
- `run_vqe_minimal_v2.py`：H₂ 最小範例（StatevectorEstimator + UCCSD，在理想模擬器求能）。
- `run_vqe_fragment_v2.py`：支援 **driver** / **integrals** 模式，並加上：
  - `--nelec`：指定總電子數或 `(nalpha,nbeta)`，以便 **Z₂ tapering**；
  - `--no-taper`：關閉 Z₂ 對稱化簡（當電子數未知時使用）。
- `make_integrals_from_driver.py`：用 **PySCF** 直接產生 **full-system** 的 `(h1,h2,enuc)`（MO 基底）。
- `make_integrals_fragment.py`：用 **PySCF+DMET** 產生 **fragment+bath 子空間** 的 `(h1,h2,enuc)`，並輸出 `frag_nelec.txt` 估算的 `nelec / (nalpha,nbeta)`。

## 安裝環境
```bash
python -m venv qiskit-env
source qiskit-env/bin/activate
pip install -r env/requirements.txt
```

## H₂ 最小範例（V2 primitives）
```bash
python3 scripts/run_vqe_minimal_v2.py
```
> V2 primitives 的 `Estimator.run([...])` 回傳值在單 observable/單參數時可能是**純量**；腳本已處理純量/陣列皆可，細節見官方文件。\[1]

## integrals 工作流
### 1) 先產生積分檔
- **full system**（MO）：
```bash
python3 scripts/make_integrals_from_driver.py   --config configs/system.yaml   --h1 frag_h1_mo.npy --h2 frag_h2_mo.npy --enuc frag_enuc.npy
```
- **fragment+bath**（DMET）：
```bash
python3 scripts/make_integrals_fragment.py   --frag-indices 0   --h1 frag_h1_mo.npy --h2 frag_h2_mo.npy --enuc frag_enuc.npy   --emit-nelec frag_nelec.txt
```

### 2) integrals 模式執行 VQE
> **重要**：若要做 **Z₂ tapering**（`problem.get_tapered_mapper`），**必須**設定 `num_particles`；否則請加 `--no-taper` 關閉。\[2]

- **指定電子數**（建議；H₂ 為封殼 2 電子）：
```bash
python3 scripts/run_vqe_fragment_v2.py --mode integrals   --h1 frag_h1_mo.npy --h2 frag_h2_mo.npy --enuc frag_enuc.npy   --nelec 2 --maxiter 200 --out results/E_fragment_v2.txt
```
- **關閉 tapering**（未提供電子數時）：
```bash
python3 scripts/run_vqe_fragment_v2.py --mode integrals   --h1 frag_h1_mo.npy --h2 frag_h2_mo.npy --enuc frag_enuc.npy --nelec 2  --no-taper --maxiter 200 --out results/E_fragment_v2.txt
```

## 參考
1. **V2 primitives**（`StatevectorEstimator`/`EstimatorV2` 的 `run([...])` 與 `pub.data.evs` 回傳格式）：IBM Quantum 官方文件。\[3]
2. **Z₂ 對稱化簡（tapering）與粒子數需求**、**Parity 映射**、`get_tapered_mapper(...)`：Qiskit Nature 映射/教學。\[4]
3. **ElectronicIntegrals** 的匯入位置（0.7 在 `second_q.operators`）與 `from_raw_integrals(h1,h2)`：Qiskit Nature 教學（Transforming Problems）。\[5]
4. **二次量子化哈密頓量不含核斥能**：需在 `ElectronicEnergy.nuclear_repulsion_energy` 另外設定；教學示例有註記。\[6]
5. **PySCF AO→MO 二體積分**：`ao2mo.full` + `ao2mo.restore(1, ...)`（還原為 4-index chemist’s notation）。\[7]\[8]

---
**註**：本 repo 採 **src-layout**；在本目錄執行腳本時，請以 repo 根目錄為工作目錄（或設定 `PYTHONPATH=src`）。

