#!/bin/bash

#DSUB -n lulcc_025deg_b__BAND_PAD__
#DSUB -A root.xibnlkjdxstbckxygcxyuan
#DSUB -q root.default
#DSUB -N 1
#DSUB -R cpu=1
#DSUB -oo /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/logs/lulcc_025deg_b__BAND_PAD___%J.out
#DSUB -eo /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/logs/lulcc_025deg_b__BAND_PAD___%J.err

set -euo pipefail

PROJECT_DIR=/home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC
CODE_DIR=${PROJECT_DIR}/Code

# Python 3.8.9
source /share/ccsuite/ENV/setenvpython389.sh
source ${PROJECT_DIR}/envs/LULCC/bin/activate

# External Python packages + project code
export PYTHONPATH="${PROJECT_DIR}/site-packages_py38:${CODE_DIR}:${PYTHONPATH:-}"

# One process / one CPU thread per band
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONUNBUFFERED=1

# 0.25 degree global bands:
# 3 degrees per band = 12 longitude cells
# 360 / 3 = 120 bands, BAND_ID = 0 ... 119
export BAND_ID=__BAND_ID__
export BAND_SIZE=3
export STAGGER_MAX=0

mkdir -p "${PROJECT_DIR}/logs"
mkdir -p "${PROJECT_DIR}/Out_ncfile"

cd "${CODE_DIR}"

echo "========================================"
echo "Start time: $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Python: $(which python)"
python --version
echo "BAND_SIZE=${BAND_SIZE}"
echo "BAND_ID=${BAND_ID}"
echo "========================================"

python -u -m server.main_025deg

echo "========================================"
echo "Finish time: $(date)"
echo "========================================"

