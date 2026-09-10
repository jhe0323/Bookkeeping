#!/bin/bash
set -euo pipefail

# Usage:
#   bash server/submit_mc_samples_1deg.sh START_SAMPLE END_SAMPLE [CPUS_PER_SAMPLE]
#
# Examples:
#   bash server/submit_mc_samples_1deg.sh 11 500
#   bash server/submit_mc_samples_1deg.sh 11 49 24

PROJECT_DIR=/home/xibnlkjdxstbckxygcxyuan/jzh0126/LULCC
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATE="${SCRIPT_DIR}/run_mc_sample_1deg.sh"

START_SAMPLE="${1:?START_SAMPLE is required}"
END_SAMPLE="${2:?END_SAMPLE is required}"
CPU_COUNT="${3:-36}"
SUBMIT_SLEEP="${SUBMIT_SLEEP:-0.05}"

if ! [[ "${START_SAMPLE}" =~ ^[0-9]+$ && "${END_SAMPLE}" =~ ^[0-9]+$ && "${CPU_COUNT}" =~ ^[0-9]+$ ]]; then
    echo "Arguments must be integers." >&2
    exit 2
fi
if (( END_SAMPLE < START_SAMPLE )); then
    echo "Invalid sample range." >&2
    exit 2
fi
if (( CPU_COUNT < 1 || CPU_COUNT > 36 )); then
    echo "CPUS_PER_SAMPLE must be 1-36." >&2
    exit 2
fi
if [[ ! -f "${TEMPLATE}" ]]; then
    echo "Template not found: ${TEMPLATE}" >&2
    exit 1
fi

GENERATED="${PROJECT_DIR}/logs/generated_mc_1deg_samples"
mkdir -p "${GENERATED}" "${PROJECT_DIR}/logs"

STAMP="$(date +%Y%m%d_%H%M%S)"
SUBMIT_LOG="${PROJECT_DIR}/logs/submit_mc_1deg_${START_SAMPLE}_${END_SAMPLE}_${STAMP}.log"

TOTAL=$((END_SAMPLE - START_SAMPLE + 1))
COUNT=0

echo "Samples: ${START_SAMPLE}-${END_SAMPLE}"
echo "Jobs: ${TOTAL}"
echo "CPUs/sample: ${CPU_COUNT}"
echo "Submission log: ${SUBMIT_LOG}"

for SAMPLE_ID in $(seq "${START_SAMPLE}" "${END_SAMPLE}"); do
    SAMPLE_PAD=$(printf "%06d" "${SAMPLE_ID}")
    JOB_SCRIPT="${GENERATED}/run_mc_1deg_s${SAMPLE_PAD}.sh"

    sed       -e "s/__SAMPLE_ID__/${SAMPLE_ID}/g"       -e "s/__SAMPLE_PAD__/${SAMPLE_PAD}/g"       -e "s/__CPU_COUNT__/${CPU_COUNT}/g"       "${TEMPLATE}" > "${JOB_SCRIPT}"

    chmod 750 "${JOB_SCRIPT}"

    echo "----- sample=${SAMPLE_ID} -----" >> "${SUBMIT_LOG}"
    if ! dsub -s "${JOB_SCRIPT}" >> "${SUBMIT_LOG}" 2>&1; then
        echo "[ERROR] dsub failed for sample=${SAMPLE_ID}; see ${SUBMIT_LOG}" >&2
        exit 1
    fi

    COUNT=$((COUNT + 1))
    if (( COUNT == 1 || COUNT % 25 == 0 || COUNT == TOTAL )); then
        echo "[SUBMITTED] ${COUNT}/${TOTAL}; latest sample=${SAMPLE_ID}"
    fi
    sleep "${SUBMIT_SLEEP}"
done

echo "[DONE] submitted ${COUNT} sample jobs."
echo "Details/job IDs: ${SUBMIT_LOG}"
echo "Use 'djob' to inspect scheduler state."
