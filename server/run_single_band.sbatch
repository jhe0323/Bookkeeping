#!/bin/bash
#SBATCH -J lulcc_band
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p intel192c
#SBATCH --nodelist=comput1
#SBATCH -o /mnt/beegfs/product/lulc0120/logs/slurm_band_%A_%a.out
#SBATCH -e /mnt/beegfs/product/lulc0120/logs/slurm_band_%A_%a.err
#SBATCH -t 72:00:00
#SBATCH --array=0-119%35

set -euo pipefail
set -x

PROJECT_DIR="/mnt/beegfs/product/lulc0120"
CODE_DIR="${PROJECT_DIR}/Code"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "${LOG_DIR}"
cd "${CODE_DIR}"

echo "Node: ${SLURM_NODELIST:-unknown}"
echo "Start: $(date)"
echo "Array job: ${SLURM_ARRAY_JOB_ID:-none} Task: ${SLURM_ARRAY_TASK_ID:-0}"

# ----------------------------------------------------------------------
# Runtime selection
#
# RESOLUTION=025deg -> server.main_025deg, default BAND_SIZE=3, 120 bands
# RESOLUTION=1deg   -> server.main_1deg, default BAND_SIZE=10, 36 bands
# ----------------------------------------------------------------------
RESOLUTION="${RESOLUTION:-025deg}"
EXPERIMENT="${EXPERIMENT:-area}"
STAGGER_MAX="${STAGGER_MAX:-600}"

case "${RESOLUTION}" in
    025deg)
        RUNNER="server.main_025deg"
        BAND_SIZE="${BAND_SIZE:-3}"
        ;;
    1deg)
        RUNNER="server.main_1deg"
        BAND_SIZE="${BAND_SIZE:-10}"
        ;;
    *)
        echo "Unsupported RESOLUTION=${RESOLUTION}; use 025deg or 1deg" >&2
        exit 2
        ;;
esac

export RESOLUTION
export EXPERIMENT
export BAND_SIZE
export BAND_ID="${SLURM_ARRAY_TASK_ID:-0}"
export STAGGER_MAX

# Prevent hidden multithreading in numerical libraries.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

# Server-local Python packages plus the Git repository.
export PYTHONPATH="${PROJECT_DIR}/site-packages:${CODE_DIR}:${PYTHONPATH:-}"

which python3
python3 -V

python3 - <<'PY'
import sys
import netCDF4
import numpy
import yaml
import pandas
import pyarrow

print("python =", sys.executable)
print("dependency check OK")
PY

LOG_FILE="${LOG_DIR}/lulcc_${RESOLUTION}_${EXPERIMENT}_${SLURM_ARRAY_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}.log"

set +e
python3 -m "${RUNNER}" \
    --experiment "${EXPERIMENT}" \
    2>&1 | tee "${LOG_FILE}"
py_status=${PIPESTATUS[0]}
set -e

echo "Python exit status: ${py_status}"

if [ "${py_status}" -ne 0 ]; then
    echo "FAILED: ${RUNNER} exited with status ${py_status}" >&2
    exit "${py_status}"
fi

echo "Done: $(date)"