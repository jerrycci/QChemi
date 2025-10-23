#!/usr/bin/env bash
# ============================================================
# Ubuntu one-click setup for MOE-equivalent open-source workflow (no conda)
# Uses Python venv + apt-based dependencies + precompiled Open Babel fix
# ============================================================

set -e
echo "=== Installing MOE-replacement quantum chemistry toolchain (venv version) ==="

# ------------------------------------------------------------
# 0️⃣  系統更新與必要依賴
# ------------------------------------------------------------
sudo apt update
sudo apt install -y build-essential cmake git wget curl python3 python3-venv python3-pip \
    libboost-all-dev libeigen3-dev libxrender1 libxext6 libsm6 libgl1-mesa-glx libglu1-mesa \
    libxml2-dev libxslt1-dev zlib1g-dev libffi-dev libssl-dev

# ------------------------------------------------------------
# 1️⃣  建立 Python venv
# ------------------------------------------------------------
VENV_DIR="$HOME/openchem_env"
if [ ! -d "$VENV_DIR" ]; then
    echo "[INFO] Creating virtual environment at $VENV_DIR"
    python3 -m venv "$VENV_DIR"
fi

source "$VENV_DIR/bin/activate"
python -m pip install --upgrade pip wheel setuptools

# ------------------------------------------------------------
# 2️⃣  安裝 RDKit (precompiled wheel)
# ------------------------------------------------------------
echo "[INFO] Installing RDKit..."
pip install "numpy<2.0" "rdkit-pypi==2022.9.5"

# ------------------------------------------------------------
# 3️⃣  安裝 PyMOL (open-source)
# ------------------------------------------------------------
echo "[INFO] Installing open-source PyMOL..."
pip install pymol-open-source

# ------------------------------------------------------------
# 4️⃣  安裝 MDAnalysis, pandas, biopython, parmed
# ------------------------------------------------------------
pip install mdanalysis pandas biopython parmed

# ------------------------------------------------------------
# 5️⃣  安裝 (APT + pip)
# ------------------------------------------------------------
echo "[INFO] Installing pdb2pqr..."
sudo apt install -y pdb2pqr

# ------------------------------------------------------------
# 6️⃣  修正版 Open Babel 安裝段落 (避免 SWIG 編譯錯誤)
# ------------------------------------------------------------
echo "[INFO] Installing Open Babel (system + wheel fix)..."
sudo apt install -y openbabel libopenbabel-dev python3-openbabel
pip install openbabel-wheel || echo "[WARN] openbabel-wheel install fallback OK"

# 驗證 Open Babel
if command -v obabel &>/dev/null; then
    echo "✅ Open Babel CLI installed: $(obabel -V)"
else
    echo "⚠️ Warning: obabel command not found — please verify Open Babel installation."
fi

# ------------------------------------------------------------
# 7️⃣  驗證版本
# ------------------------------------------------------------
echo "=== Checking tool versions ==="
python -c "import rdkit; print('RDKit:', rdkit.__version__)"
pdb2pqr --version || echo "pdb2pqr OK"
python -c "import openbabel; print('OpenBabel Python API OK')" || echo "[WARN] openbabel import test skipped"

echo
echo "✅ Installation complete!"
echo "Activate environment with:"
echo "  source ~/openchem_env/bin/activate"
