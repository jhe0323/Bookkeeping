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
CODE_DIR=${PROJECT_DIR}/Code

INPUT_DIR=${PROJECT_DIR}/Out_ncfile/1deg/baseline_1deg_voidpft_v1
FINAL_OUTPUT=${INPUT_DIR}/summary_1deg.global_time_ge1850.nc

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

NFILES=$(find "${INPUT_DIR}" -maxdepth 1 -name 'summary_1deg.rank*.nc' | wc -l)
echo "Band files found: ${NFILES}"
if [[ "${NFILES}" -ne 36 ]]; then
    echo "ERROR: expected 36 files, found ${NFILES}" >&2
    exit 1
fi

TMP_WORK=$(mktemp -d /tmp/lulcc_merge_${USER}_XXXXXX)
echo "Temporary directory: ${TMP_WORK}"
cleanup() { rm -rf "${TMP_WORK}"; }
trap cleanup EXIT

# Copy each rank only once from shared storage; the Python merger then verifies
# done==1, coordinates, schemas and run/model metadata before writing output.
echo "Copying rank files to local temporary disk..."
cp "${INPUT_DIR}"/summary_1deg.rank*.nc "${TMP_WORK}/"

LOCAL_OUTPUT=${TMP_WORK}/summary_1deg.global_time_ge1850.nc
python -u "${CODE_DIR}/tools/merge_bands.py" \
    --pattern "${TMP_WORK}/summary_1deg.rank*.nc" \
    --time-min 1850 \
    --time-block 32 \
    --compression 1 \
    --output "${LOCAL_OUTPUT}" \
    --overwrite

cp "${LOCAL_OUTPUT}" "${FINAL_OUTPUT}"
ls -lh "${FINAL_OUTPUT}"
echo "Finished: $(date)"
