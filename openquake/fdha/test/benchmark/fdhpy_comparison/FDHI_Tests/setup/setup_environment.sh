#!/bin/bash
# =============================================================================
# Setup script for PFDHA comparison tests
# This script sets up the environment for running comparison tests between
# fdhpy (reference) and pfdha (this implementation)
# =============================================================================

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
COMPARISON_DIR="$PROJECT_ROOT/comparison_tests"

echo "=============================================="
echo "PFDHA Comparison Tests - Environment Setup"
echo "=============================================="
echo ""

# Check Python version
echo "[1/5] Checking Python version..."
python_version=$(python3 --version 2>&1 | awk '{print $2}')
echo "      Python version: $python_version"

# Create virtual environment (optional)
read -p "Do you want to create a new virtual environment? [y/N] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "[2/5] Creating virtual environment..."
    python3 -m venv "$COMPARISON_DIR/venv"
    source "$COMPARISON_DIR/venv/bin/activate"
    pip install --upgrade pip
else
    echo "[2/5] Skipping virtual environment creation..."
fi

# Install requirements
echo "[3/5] Installing requirements..."
pip install -r "$SCRIPT_DIR/requirements.txt"

# Install pfdha in development mode
echo "[4/5] Installing pfdha in development mode..."
pip install -e "$PROJECT_ROOT"

# Check if fdhpy is available
echo "[5/5] Checking fdhpy installation..."
if python3 -c "import fdhpy; print(f'fdhpy version: {fdhpy.__version__}')" 2>/dev/null; then
    echo "      fdhpy is installed correctly."
else
    echo "      WARNING: fdhpy is not installed."
    echo "      Please install it with: pip install fdhpy"
    echo "      Or from local clone: pip install -e /path/to/fdhpy"
fi

echo ""
echo "=============================================="
echo "Setup complete!"
echo ""
echo "To run the comparison tests:"
echo "  cd $COMPARISON_DIR"
echo "  pytest tests/ -v --html=reports/test_report.html"
echo "=============================================="

