#!/bin/bash

#DSUB -n lulcc_1deg_b000
#DSUB -A root.xibnlkjdxstbckxygcxyuan
#DSUB -q root.default
#DSUB -N 1
#DSUB -R cpu=1
#DSUB -oo /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/logs/lulcc_1deg_b000_%J.out
#DSUB -eo /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/logs/lulcc_1deg_b000_%J.err

set -euo pipefail

PROJECT_DIR=/home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC
CODE_DIR=${PROJECT_DIR}/Code

# 加载平台 Python 3.8.9
source /share/ccsuite/ENV/setenvpython389.sh

# 激活模型虚拟环境
source ${PROJECT_DIR}/envs/LULCC/bin/activate

# 避免数值库在单核作业中额外启动大量线程
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1

# 模型运行参数
export RESOLUTION=1deg
export EXPERIMENT=area
export BAND_SIZE=10
export BAND_ID=0
export STAGGER_MAX=0

export HDF5_USE_FILE_LOCKING=FALSE
# 当前只测试第 0 个 10° 经度条带
export BAND_ID=0
export BAND_SIZE=10
export STAGGER_MAX=0
export PYTHONPATH="${PROJECT_DIR}/site-packages_py38:${CODE_DIR}:${PYTHONPATH:-}"

mkdir -p "${PROJECT_DIR}/logs"
mkdir -p "${PROJECT_DIR}/Out_ncfile"

cd "${CODE_DIR}"

echo "========================================"
echo "Job ID: ${CCSCHEDULER_JOB_ID:-unknown}"
echo "Start time: $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Python: $(which python)"
python --version
#echo "RESOLUTION=${RESOLUTION}"
#echo "EXPERIMENT=${EXPERIMENT}"
echo "BAND_SIZE=${BAND_SIZE}"
echo "BAND_ID=${BAND_ID}"
echo "========================================"

python -u -m server.main_1deg

echo "========================================"
echo "Finish time: $(date)"
echo "========================================"
