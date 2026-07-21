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

PROJECT_DIR="/mnt/beegfs/product/lulc0120"
CODE_DIR="${PROJECT_DIR}/Code"
RESOLUTION="${RESOLUTION:-025deg}"

case "${RESOLUTION}" in
    025deg) RUNNER="server.main_025deg" ;;
    1deg)   RUNNER="server.main_1deg" ;;
    *) echo "Unsupported RESOLUTION=${RESOLUTION}" >&2; exit 2 ;;
esac

mkdir -p "${PROJECT_DIR}/logs"
cd "${CODE_DIR}"
export BAND_ID="${SLURM_ARRAY_TASK_ID:-0}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONPATH="${PROJECT_DIR}/site-packages:${CODE_DIR}:${PYTHONPATH:-}"

python3 -m "${RUNNER}"
