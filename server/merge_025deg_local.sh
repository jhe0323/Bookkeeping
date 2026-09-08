#!/bin/bash

#DSUB -n merge_025deg
#DSUB -A root.xibnlkjdxstbckxygcxyuan
#DSUB -q root.default
#DSUB -N 1
#DSUB -R cpu=1
#DSUB -oo /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/logs/merge_025deg_%J.out
#DSUB -eo /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/logs/merge_025deg_%J.err

set -euo pipefail

PROJECT_DIR=/home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC
CODE_DIR=${PROJECT_DIR}/Code

INPUT_DIR=${PROJECT_DIR}/Out_ncfile/025deg/baseline_025deg_v1
FINAL_OUTPUT=${INPUT_DIR}/summary_025deg.global_time_ge1850.nc

source /share/ccsuite/ENV/setenvpython389.sh
source ${PROJECT_DIR}/envs/LULCC/bin/activate

export PYTHONPATH="${PROJECT_DIR}/site-packages_py38:${CODE_DIR}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONUNBUFFERED=1

mkdir -p "${PROJECT_DIR}/logs"

NFILES=$(find "${INPUT_DIR}" -maxdepth 1 -name 'summary_025deg.rank*.nc' | wc -l)
echo "Band files found: ${NFILES}"
if [[ "${NFILES}" -ne 120 ]]; then
    echo "ERROR: expected 120 files, found ${NFILES}" >&2
    exit 1
fi

# 0.25-degree rank files are merged directly on shared storage. The Python
# tool streams time blocks and verifies every rank before producing the output.
python -u "${CODE_DIR}/tools/merge_bands.py" \
    --config "${CODE_DIR}/config/run_025deg.yml" \
    --pattern "${INPUT_DIR}/summary_025deg.rank*.nc" \
    --time-min 1850 \
    --time-block 32 \
    --compression 1 \
    --output "${FINAL_OUTPUT}" \
    --overwrite

ls -lh "${FINAL_OUTPUT}"
echo "Finished: $(date)"
