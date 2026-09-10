#!/bin/bash
# One DSUB job = one MC sample; all 36 bands run inside this job.

#DSUB -n mc1d_s__SAMPLE_PAD__
#DSUB -A root.xibnlkjdxstbckxygcxyuan
#DSUB -q root.default
#DSUB -N 1
#DSUB -R cpu=__CPU_COUNT__
#DSUB -oo /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/logs/mc1d_s__SAMPLE_PAD___%J.out
#DSUB -eo /home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC/logs/mc1d_s__SAMPLE_PAD___%J.err

set -uo pipefail

PROJECT_DIR=/home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC
CODE_DIR=${PROJECT_DIR}/Code
source "${PROJECT_DIR}/source.txt"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONUNBUFFERED=1

export MC_SAMPLE_ID=__SAMPLE_ID__
export BAND_SIZE=10
export STAGGER_MAX=0
export RUN_CONFIG="${CODE_DIR}/config/run_mc_1deg.yml"

MAX_PARALLEL=__CPU_COUNT__
SAMPLE_PAD=__SAMPLE_PAD__
RUN_TAG="$(date +%Y%m%d_%H%M%S)_$$"
BAND_LOG_DIR="${PROJECT_DIR}/logs/mc_1deg_band_logs/sample_${SAMPLE_PAD}/${RUN_TAG}"

mkdir -p "${BAND_LOG_DIR}" "${PROJECT_DIR}/Out_ncfile/MC"
cd "${CODE_DIR}"

echo "Start: $(date)"
echo "Host: $(hostname)"
echo "Sample: ${MC_SAMPLE_ID}"
echo "Bands: 0-35"
echo "Parallel workers: ${MAX_PARALLEL}"
echo "Band logs: ${BAND_LOG_DIR}"

FAILED_BANDS=()

wait_wave() {
    local item pid band status
    for item in "$@"; do
        pid="${item%%:*}"
        band="${item##*:}"
        if wait "${pid}"; then
            :
        else
            status=$?
            echo "[FAIL] sample=${MC_SAMPLE_ID} band=${band} exit=${status}"
            FAILED_BANDS+=("${band}")
        fi
    done
}

WAVE=()

for BAND_ID in $(seq 0 35); do
    BAND_PAD=$(printf "%03d" "${BAND_ID}")
    (
        export BAND_ID="${BAND_ID}"
        python -u -m server.main_mc_1deg
    ) >"${BAND_LOG_DIR}/band_${BAND_PAD}.out" 2>"${BAND_LOG_DIR}/band_${BAND_PAD}.err" &

    WAVE+=("$!:${BAND_ID}")

    if (( ${#WAVE[@]} >= MAX_PARALLEL )); then
        wait_wave "${WAVE[@]}"
        WAVE=()
    fi
done

if (( ${#WAVE[@]} > 0 )); then
    wait_wave "${WAVE[@]}"
fi

if (( ${#FAILED_BANDS[@]} > 0 )); then
    echo "[FAILED] sample=${MC_SAMPLE_ID}; bands=${FAILED_BANDS[*]}"
    echo "Inspect ${BAND_LOG_DIR}"
    exit 1
fi

echo "[PASS] sample=${MC_SAMPLE_ID}; all bands completed or were skipped as current."
echo "Finish: $(date)"
