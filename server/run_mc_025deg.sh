#!/bin/bash

#DSUB -n mc_025deg_s__SAMPLE_PAD___b__BAND_PAD__
#DSUB -A root.xibnlkjdxstbckxygcxyuan
#DSUB -q root.default
#DSUB -N 1
#DSUB -R cpu=1
#DSUB -oo /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/logs/mc_025deg_s__SAMPLE_PAD___b__BAND_PAD___%J.out
#DSUB -eo /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/logs/mc_025deg_s__SAMPLE_PAD___b__BAND_PAD___%J.err

set -euo pipefail

PROJECT_DIR=/home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC
CODE_DIR=${PROJECT_DIR}/Code

source /share/ccsuite/ENV/setenvpython389.sh
source ${PROJECT_DIR}/envs/LULCC/bin/activate

export PYTHONPATH="${PROJECT_DIR}/site-packages_py38:${CODE_DIR}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONUNBUFFERED=1

export MC_SAMPLE_ID=__SAMPLE_ID__
export BAND_ID=__BAND_ID__
export BAND_SIZE=3
export STAGGER_MAX=0
export RUN_CONFIG="${CODE_DIR}/config/run_mc_025deg.yml"

mkdir -p "${PROJECT_DIR}/logs"
mkdir -p "${PROJECT_DIR}/Out_ncfile/MC"

cd "${CODE_DIR}"

echo "========================================"
echo "Start time: $(date)"
echo "Host: $(hostname)"
echo "Python: $(which python)"
python --version
echo "Resolution: 025deg"
echo "MC_SAMPLE_ID=${MC_SAMPLE_ID}"
echo "BAND_ID=${BAND_ID} / 119"
echo "RUN_CONFIG=${RUN_CONFIG}"
echo "========================================"

python -u -m server.main_mc_025deg

echo "========================================"
echo "Finish time: $(date)"
echo "========================================"
