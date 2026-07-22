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

# The array range is intentionally NOT defined in this file.
# It must be supplied explicitly when submitting, for example:
#
# 0.25 degree:
#   sbatch --array=0-119%35 --export=ALL,RESOLUTION=025deg server/run_single_band.sbatch
#
# 1 degree:
#   sbatch --array=0-35%25 --export=ALL,RESOLUTION=1deg server/run_single_band.sbatch
#
# Single-band rerun:
#   sbatch --array=31-31 --export=ALL,RESOLUTION=025deg server/run_single_band.sbatch

set -euo pipefail

PROJECT_DIR="/mnt/beegfs/product/lulc0120"
CODE_DIR="${PROJECT_DIR}/Code"
RESOLUTION="${RESOLUTION:-025deg}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "This script must be submitted with an explicit --array range." >&2
    echo "Examples:" >&2
    echo "  025deg: sbatch --array=0-119%35 --export=ALL,RESOLUTION=025deg server/run_single_band.sbatch" >&2
    echo "  1deg:   sbatch --array=0-35%25  --export=ALL,RESOLUTION=1deg   server/run_single_band.sbatch" >&2
    exit 2
fi

case "${RESOLUTION}" in
    025deg) RUNNER="server.main_025deg" ;;
    1deg)   RUNNER="server.main_1deg" ;;
    *) echo "Unsupported RESOLUTION=${RESOLUTION}" >&2; exit 2 ;;
esac

mkdir -p "${PROJECT_DIR}/logs"
cd "${CODE_DIR}"
export BAND_ID="${SLURM_ARRAY_TASK_ID}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONPATH="${PROJECT_DIR}/site-packages:${CODE_DIR}:${PYTHONPATH:-}"

python3 -m "${RUNNER}"
