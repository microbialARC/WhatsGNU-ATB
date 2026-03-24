#!/usr/bin/env bash
set -euo pipefail

# ── WhatsGNU-ATB Setup Script ──
# Creates conda environment, installs dependencies, and clones the repo.
#
# Usage:
#   bash setup_whatsgnu_atb.sh
#   bash setup_whatsgnu_atb.sh --env-name my_env --install-dir /path/to/install

ENV_NAME="whatsgnu-atb"
INSTALL_DIR="$(pwd)/WhatsGNU-ATB"
PYTHON_VERSION="3.12"
REPO_URL="https://github.com/microbialARC/WhatsGNU-ATB.git"

# ── Parse arguments ──
while [[ $# -gt 0 ]]; do
    case "$1" in
        --env-name)   ENV_NAME="$2"; shift 2 ;;
        --install-dir) INSTALL_DIR="$2"; shift 2 ;;
        --python)     PYTHON_VERSION="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: bash setup_whatsgnu_atb.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --env-name NAME    Conda environment name (default: whatsgnu-atb)"
            echo "  --install-dir DIR  Installation directory (default: ./WhatsGNU-ATB)"
            echo "  --python VERSION   Python version (default: 3.12)"
            echo "  -h, --help         Show this help"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

echo "============================================"
echo "  WhatsGNU-ATB Setup"
echo "============================================"
echo "  Conda env:    $ENV_NAME"
echo "  Install dir:  $INSTALL_DIR"
echo "  Python:       $PYTHON_VERSION"
echo "  Repo:         $REPO_URL"
echo ""

# ── Check conda ──
if ! command -v conda &>/dev/null; then
    echo "ERROR: conda not found. Install Miniconda/Anaconda first."
    exit 1
fi

# ── Create conda environment ──
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "Conda environment '$ENV_NAME' already exists. Skipping creation."
else
    echo "Creating conda environment '$ENV_NAME' with Python $PYTHON_VERSION ..."
    conda create -n "$ENV_NAME" -c conda-forge python="$PYTHON_VERSION" -y
fi

echo "Installing dependencies ..."
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"
pip install numpy lmdb pandas matplotlib seaborn networkx adjustText scipy

echo ""
echo "Installed versions:"
python --version
python -c "import numpy; print(f'  numpy {numpy.__version__}')"
python -c "import lmdb; print(f'  lmdb {lmdb.version()}')"
python -c "import pandas; print(f'  pandas {pandas.__version__}')"

# ── Clone repo ──
if [ -d "$INSTALL_DIR" ]; then
    echo ""
    echo "Directory $INSTALL_DIR already exists. Pulling latest ..."
    cd "$INSTALL_DIR"
    git pull
else
    echo ""
    echo "Cloning WhatsGNU-ATB ..."
    git clone "$REPO_URL" "$INSTALL_DIR"
    cd "$INSTALL_DIR"
fi

chmod +x scripts/*.py 2>/dev/null || true

echo ""
echo "============================================"
echo "  Setup complete!"
echo "============================================"
echo ""
echo "To activate the environment:"
echo "  conda activate $ENV_NAME"
echo ""
echo "To download the pre-built database:"
echo "  python $INSTALL_DIR/scripts/download_osf.py --folder WGNU_ATB_DB --out-dir ./WGNU_ATB_DB"
echo ""
echo "To query a genome:"
echo "  python $INSTALL_DIR/scripts/Query_WhatsGNU_ATB.py \\"
echo "      --db_dir ./WGNU_ATB_DB \\"
echo "      --shards 8 \\"
echo "      --faa your_genome.bakta.faa \\"
echo "      --out_dir results/"
