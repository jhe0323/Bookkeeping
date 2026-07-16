#!/bin/bash
set -euo pipefail

# Usage:
#   bash server/clean_run.sh 025deg area
#   bash server/clean_run.sh 1deg bio-strict
#   bash server/clean_run.sh all
#
# This script never touches In_ncfile/, archive/, site-packages/,
# outputs/, or outputs_1deg/.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PROJECT_DIR="$(cd "${CODE_DIR}/.." && pwd)"

OUT_ROOT="${PROJECT_DIR}/Out_ncfile"
LOG_DIR="${PROJECT_DIR}/logs"

TARGET="${1:-}"
EXPERIMENT_ALIAS="${2:-area}"

usage() {
    echo "Usage:"
    echo "  bash server/clean_run.sh 025deg {area|bio-strict|bio-forced}"
    echo "  bash server/clean_run.sh 1deg   {area|bio-strict|bio-forced}"
    echo "  bash server/clean_run.sh all"
}

case "${EXPERIMENT_ALIAS}" in
    area)
        EXPERIMENT_NAME="harvest_area_driven"
        ;;
    bio-strict)
        EXPERIMENT_NAME="harvest_bio_strict"
        ;;
    bio-forced)
        EXPERIMENT_NAME="harvest_bio_forced"
        ;;
    *)
        echo "Unknown experiment alias: ${EXPERIMENT_ALIAS}" >&2
        usage
        exit 2
        ;;
esac

case "${TARGET}" in
    025deg|1deg)
        OUTPUT_TARGET="${OUT_ROOT}/${TARGET}/${EXPERIMENT_NAME}"
        echo "The following will be removed:"
        echo "  output files: ${OUTPUT_TARGET}/summary_*.nc"
        echo "  selected logs: ${LOG_DIR}/lulcc_${TARGET}_${EXPERIMENT_ALIAS}_*.log"
        echo "  SLURM logs: ${LOG_DIR}/slurm_band_*.out and *.err"
        ;;
    all)
        OUTPUT_TARGET="${OUT_ROOT}"
        echo "The following will be removed:"
        echo "  all files below: ${OUT_ROOT}"
        echo "  all files below: ${LOG_DIR}"
        ;;
    *)
        usage
        exit 2
        ;;
esac

read -r -p "Confirm deletion? (y/N): " confirm
if [[ "${confirm}" != "y" && "${confirm}" != "Y" ]]; then
    echo "Cancelled."
    exit 0
fi

if [[ "${TARGET}" == "all" ]]; then
    if [[ -d "${OUT_ROOT}" ]]; then
        find "${OUT_ROOT}" -mindepth 1 -type f -delete
        find "${OUT_ROOT}" -mindepth 1 -type d -empty -delete
    fi
    if [[ -d "${LOG_DIR}" ]]; then
        find "${LOG_DIR}" -mindepth 1 -type f -delete
    fi
else
    if [[ -d "${OUTPUT_TARGET}" ]]; then
        find "${OUTPUT_TARGET}" -maxdepth 1 -type f -name 'summary_*.nc' -delete
    fi

    rm -f \
        "${LOG_DIR}/lulcc_${TARGET}_${EXPERIMENT_ALIAS}_"*.log \
        "${LOG_DIR}/slurm_band_"*.out \
        "${LOG_DIR}/slurm_band_"*.err
fi

echo "Cleanup complete."
