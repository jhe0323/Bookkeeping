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
TOOLS_DIR=${PROJECT_DIR}/Tools

INPUT_DIR=${PROJECT_DIR}/Out_ncfile/025deg/baseline_025deg_v1
FINAL_OUTPUT=${INPUT_DIR}/summary_025deg.global_time_ge1850.nc

source /share/ccsuite/ENV/setenvpython389.sh
source ${PROJECT_DIR}/envs/LULCC/bin/activate

export PYTHONPATH="${PROJECT_DIR}/site-packages_py38:${PROJECT_DIR}/Code:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONUNBUFFERED=1

echo "========================================"
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "Python: $(which python)"
python --version
echo "Input: ${INPUT_DIR}"
echo "Output: ${FINAL_OUTPUT}"
echo "========================================"

# 检查120个rank
NFILES=$(find "${INPUT_DIR}" -maxdepth 1 \
    -name 'summary_025deg.rank*.nc' | wc -l)

echo "Band files found: ${NFILES}"

if [[ "${NFILES}" -ne 120 ]]; then
    echo "ERROR: expected 120 files, found ${NFILES}" >&2
    exit 1
fi

echo "Input size:"
du -sh "${INPUT_DIR}"

echo "Home filesystem:"
df -h "${INPUT_DIR}"

echo
echo "Starting merge directly on home storage..."
date

python -u "${TOOLS_DIR}/2.merge_bands.py" \
    --pattern "${INPUT_DIR}/summary_025deg.rank*.nc" \
    --time-min 1850 \
    --time-block 64 \
    --compression 1 \
    --output "${FINAL_OUTPUT}" \
    --overwrite

echo
echo "Merge finished."
date

ls -lh "${FINAL_OUTPUT}"

echo "========================================"
echo "Finished: $(date)"
echo "========================================"
