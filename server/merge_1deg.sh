#!/bin/bash

#DSUB -n merge_1deg
#DSUB -A root.xibnlkjdxstbckxygcxyuan
#DSUB -q root.default
#DSUB -N 1
#DSUB -R cpu=1
#DSUB -oo /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/logs/merge_1deg_%J.out
#DSUB -eo /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/logs/merge_1deg_%J.err

set -euo pipefail

PROJECT_DIR=/home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC
TOOLS_DIR=${PROJECT_DIR}/Tools

INPUT_DIR=${PROJECT_DIR}/Out_ncfile/1deg/baseline_1deg_voidpft_v1/
FINAL_OUTPUT=${INPUT_DIR}/summary_1deg.global_time_ge1850.nc

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
echo "========================================"

# --------------------------------------------------
# 1. 检查36个输入文件
# --------------------------------------------------

NFILES=$(find "${INPUT_DIR}" -maxdepth 1 \
    -name 'summary_1deg.rank*.nc' | wc -l)

echo "Band files found: ${NFILES}"

if [[ "${NFILES}" -ne 36 ]]; then
    echo "ERROR: expected 36 files, found ${NFILES}" >&2
    exit 1
fi

# --------------------------------------------------
# 2. 创建计算节点本地临时目录
# --------------------------------------------------

TMP_WORK=$(mktemp -d /tmp/lulcc_merge_${USER}_XXXXXX)

echo "Temporary directory:"
echo "  ${TMP_WORK}"

df -h /tmp

cleanup() {
    echo "Cleaning temporary files..."
    rm -rf "${TMP_WORK}"
}
trap cleanup EXIT

# --------------------------------------------------
# 3. rank文件只从共享存储复制一次
# --------------------------------------------------

echo "Copying rank files to local temporary disk..."
date

cp "${INPUT_DIR}"/summary_1deg.rank*.nc "${TMP_WORK}/"

echo "Copy finished."
date

du -sh "${TMP_WORK}"

# --------------------------------------------------
# 4. 在本地磁盘完成merge
# --------------------------------------------------

LOCAL_OUTPUT=${TMP_WORK}/summary_1deg.global_time_ge1850.nc

echo "Starting local merge..."
date

python -u "${TOOLS_DIR}/2.merge_bands.py" \
    --pattern "${TMP_WORK}/summary_1deg.rank*.nc" \
    --time-min 1850 \
    --time-block 32 \
    --compression 1 \
    --output "${LOCAL_OUTPUT}" \
    --overwrite

echo "Local merge finished."
date

ls -lh "${LOCAL_OUTPUT}"

# --------------------------------------------------
# 5. 最终结果只写回共享存储一次
# --------------------------------------------------

echo "Copying final output back..."
date

cp "${LOCAL_OUTPUT}" "${FINAL_OUTPUT}"

echo "Copy finished."
date

ls -lh "${FINAL_OUTPUT}"

echo "========================================"
echo "Finished: $(date)"
echo "========================================"
