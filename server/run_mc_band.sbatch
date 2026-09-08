#!/bin/bash
#SBATCH -J lulcc_mc
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p intel192c
#SBATCH --nodelist=comput1
#SBATCH -o /mnt/beegfs/product/lulc0120/logs/mc_%A_%a.out
#SBATCH -e /mnt/beegfs/product/lulc0120/logs/mc_%A_%a.err
#SBATCH -t 72:00:00

# This job maps one SLURM array element to one (sample_id, band_id) pair.
# Submit through server/submit_mc_batch.sh rather than specifying --array by hand.

set -euo pipefail

PROJECT_DIR="/mnt/beegfs/product/lulc0120"
CODE_DIR="${PROJECT_DIR}/Code"
RUN_CONFIG="${RUN_CONFIG:-${CODE_DIR}/config/run_mc_1deg.yml}"
BANDS_TOTAL="${BANDS_TOTAL:-36}"
MC_SAMPLE_OFFSET="${MC_SAMPLE_OFFSET:-0}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "This script must run as a SLURM array job." >&2
    exit 2
fi

LOCAL_SAMPLE=$(( SLURM_ARRAY_TASK_ID / BANDS_TOTAL ))
BAND_ID=$(( SLURM_ARRAY_TASK_ID % BANDS_TOTAL ))
MC_SAMPLE_ID=$(( MC_SAMPLE_OFFSET + LOCAL_SAMPLE ))

mkdir -p "${PROJECT_DIR}/logs"
cd "${CODE_DIR}"

export MC_SAMPLE_ID
export BAND_ID
export RUN_CONFIG
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONPATH="${PROJECT_DIR}/site-packages:${CODE_DIR}:${PYTHONPATH:-}"

echo "[MC JOB] task=${SLURM_ARRAY_TASK_ID} sample=${MC_SAMPLE_ID} band=${BAND_ID}/${BANDS_TOTAL}"
python3 -m server.main_mc_1deg
