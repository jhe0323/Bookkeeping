#!/bin/bash
#SBATCH -J lulcc_band
#SBATCH -N 1
#SBATCH --ntasks=1                 # 单任务
#SBATCH --cpus-per-task=1          # 单核
#SBATCH -p intel192c
#SBATCH --nodelist=comput1         # 可按需注释
#SBATCH -o slurm_band_%A_%a.out
#SBATCH -e slurm_band_%A_%a.err
#SBATCH -t 72:00:00
#SBATCH --array=0-119%35

set -x
set -o pipefail
cd "$SLURM_SUBMIT_DIR"

echo "Node: $SLURM_NODELIST"
echo "Start: $(date)"
echo "Array job: $SLURM_ARRAY_JOB_ID  Task: $SLURM_ARRAY_TASK_ID"

# ========= 1) 统一环境 =========
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

# Python 包路径（按你的环境）
export PYTHONPATH="/mnt/beegfs/product/lulc0120/site-packages:${PYTHONPATH:-}"
# ========= 2) 分带参数 =========
export BAND_SIZE=3 #修改带宽，并修改相应LN11 纬度范围
export BAND_ID=$SLURM_ARRAY_TASK_ID

# ========= 3) 依赖自检（可选） =========
which python3
python3 -V
python3 - <<'PY' || exit 1
import sys, site
print("python =", sys.executable)
print("usersite =", getattr(site, "getusersitepackages", lambda: "N/A")())
import netCDF4, numpy, yaml
print("netCDF4/numpy/yaml OK")
import pandas
print("pandas OK")

import pyarrow
print("pyarrow OK")
PY

# ========= 4) 运行（无需 mpirun） =========
mkdir -p logs
set +e
python3 main.py 2>&1 | tee "logs/lulcc_band_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
py_status=${PIPESTATUS[0]}
set -e

echo "Python exit status: ${py_status}"

if [ "${py_status}" -ne 0 ]; then
    echo "FAILED: python3 main.py exited with status ${py_status}"
    exit "${py_status}"
fi

echo "Done: $(date)"
