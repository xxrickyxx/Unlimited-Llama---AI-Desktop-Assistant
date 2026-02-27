#!/usr/bin/env bash
# run_experiment.sh — Run the Goldbach fiber-sum experiment with no further user input.
#
# Usage (all arguments are optional; defaults are N=1000000, W=1024, A=2):
#   bash run_experiment.sh [--N N] [--W W] [--A A] [--save-csv PATH] [--quiet]
#
# Examples:
#   bash run_experiment.sh
#   bash run_experiment.sh --A 3
#   bash run_experiment.sh --N 50000 --W 256 --A 2 --save-csv results.csv

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Ensure numpy is available.
if ! python3 -c "import numpy" 2>/dev/null; then
    echo "[run_experiment] Installing numpy..."
    # Use --user when not inside a virtual/conda environment to avoid
    # modifying the system-wide Python installation.
    if [ -n "${VIRTUAL_ENV:-}" ] || [ -n "${CONDA_PREFIX:-}" ]; then
        pip install --quiet numpy
    else
        pip install --quiet --user numpy
    fi
fi

# Default to A=2 so minor-arc sampling works out of the box at N=1_000_000.
# Only inject the default when the caller has not already supplied --A.
ARGS=()
if [[ ! " $* " =~ " --A " ]]; then
    ARGS+=("--A" "2")
fi
ARGS+=("$@")

echo "[run_experiment] Starting Goldbach fiber-sum experiment (adesso!)..."
python3 "${SCRIPT_DIR}/goldbach_experiment.py" "${ARGS[@]}"
