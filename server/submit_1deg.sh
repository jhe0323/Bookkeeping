#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TEMPLATE="${SCRIPT_DIR}/run_1deg.sh"
GENERATED="${SCRIPT_DIR}/generated_1deg"

START_BAND="${1:-0}"
END_BAND="${2:-35}"

if [[ ! -f "${TEMPLATE}" ]]; then
    echo "Template not found:"
    echo "  ${TEMPLATE}"
    exit 1
fi

if (( START_BAND < 0 || END_BAND > 35 || START_BAND > END_BAND )); then
    echo "Invalid band range: ${START_BAND}-${END_BAND}"
    echo "Allowed range: 0-35"
    exit 2
fi

mkdir -p "${GENERATED}"

echo "========================================"
echo "Submitting 1deg LULCC bands"
echo "Range: ${START_BAND}-${END_BAND}"
echo "Generated scripts:"
echo "  ${GENERATED}"
echo "========================================"

for BAND_ID in $(seq "${START_BAND}" "${END_BAND}"); do

    BAND_PAD=$(printf "%03d" "${BAND_ID}")

    JOB_SCRIPT="${GENERATED}/run_1deg_b${BAND_PAD}.sh"

    sed \
        -e "s/__BAND_ID__/${BAND_ID}/g" \
        -e "s/__BAND_PAD__/${BAND_PAD}/g" \
        "${TEMPLATE}" > "${JOB_SCRIPT}"

    chmod 750 "${JOB_SCRIPT}"

    echo
    echo "Submitting BAND_ID=${BAND_ID} ..."

    dsub -s "${JOB_SCRIPT}"

    # 避免瞬间大量请求调度器
    sleep 0.2
done

echo
echo "========================================"
echo "Submission complete."
echo "Bands: ${START_BAND}-${END_BAND}"
echo "========================================"
