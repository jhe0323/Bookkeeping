#!/bin/bash

#DSUB -n summarize_025deg
#DSUB -A root.xibnlkjdxstbckxygcxyuan
#DSUB -q root.default
#DSUB -N 1
#DSUB -R cpu=1
#DSUB -oo /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/logs/summarize_025deg_%J.out
#DSUB -eo /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/logs/summarize_025deg_%J.err

set -euo pipefail

PROJECT_DIR=/home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC
TOOLS_DIR=${PROJECT_DIR}/Tools
INPUT_DIR=${PROJECT_DIR}/Out_ncfile/025deg/baseline_025deg_v1
TIME_MIN=${TIME_MIN:-1850}
OUTPUT=${INPUT_DIR}/summary_025deg.global_time_ge${TIME_MIN}_dig.csv

source /share/ccsuite/ENV/setenvpython389.sh
source ${PROJECT_DIR}/envs/LULCC/bin/activate
export PYTHONPATH="${PROJECT_DIR}/site-packages_py38:${PROJECT_DIR}/Code:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONUNBUFFERED=1

echo "Start: $(date)"
echo "Host: $(hostname)"
echo "Input: ${INPUT_DIR}"
echo "Time min: ${TIME_MIN}"
echo "Output: ${OUTPUT}"

python -u "${TOOLS_DIR}/5.NC_bands_to_csv.py" \
    --pattern "${INPUT_DIR}/summary_025deg.rank*.nc" \
    --output "${OUTPUT}" \
    --time-min "${TIME_MIN}" \
    --expected-bands 120 \
    --progress-every 1

echo "Finished: $(date)"
ls -lh "${OUTPUT}"
