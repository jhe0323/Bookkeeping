#!/bin/bash
set -euo pipefail

# Clean model outputs generated under:
#   Out_ncfile/<resolution>/<run.name>/
#
# Usage:
#   bash server/clean_run.sh 025deg baseline_025deg_area
#   bash server/clean_run.sh 1deg baseline_1deg_area
#   bash server/clean_run.sh 1deg local_1deg_test
#   bash server/clean_run.sh all
#
# Optional:
#   CLEAN_LOGS=1 bash server/clean_run.sh 025deg baseline_025deg_area
#
# CLEAN_LOGS=1 removes SLURM output/error logs after deleting the selected run.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PROJECT_DIR="$(cd "${CODE_DIR}/.." && pwd)"

OUT_ROOT="${PROJECT_DIR}/Out_ncfile"
LOG_DIR="${PROJECT_DIR}/logs"

TARGET="${1:-}"
RUN_NAME="${2:-}"
CLEAN_LOGS="${CLEAN_LOGS:-0}"


usage() {
    echo "Usage:"
    echo "  bash server/clean_run.sh 025deg <run.name>"
    echo "  bash server/clean_run.sh 1deg   <run.name>"
    echo "  bash server/clean_run.sh all"
    echo
    echo "Examples:"
    echo "  bash server/clean_run.sh 025deg baseline_025deg_area"
    echo "  bash server/clean_run.sh 1deg baseline_1deg_area"
    echo "  CLEAN_LOGS=1 bash server/clean_run.sh 1deg baseline_1deg_area"
}


validate_run_name() {
    local name="$1"

    if [[ -z "${name}" ]]; then
        echo "Error: run.name is required." >&2
        usage
        exit 2
    fi

    # Prevent paths such as ../, /tmp/test, or nested directories.
    if [[ "${name}" == "." || "${name}" == ".." || "${name}" == *"/"* ]]; then
        echo "Error: invalid run.name: ${name}" >&2
        exit 2
    fi
}


remove_logs() {
    if [[ "${CLEAN_LOGS}" != "1" ]]; then
        return
    fi

    if [[ ! -d "${LOG_DIR}" ]]; then
        echo "Log directory does not exist: ${LOG_DIR}"
        return
    fi

    echo "Removing SLURM logs from:"
    echo "  ${LOG_DIR}"

    find "${LOG_DIR}" -maxdepth 1 -type f \
        \( \
            -name 'slurm_band_*.out' \
            -o -name 'slurm_band_*.err' \
            -o -name 'slurm-*.out' \
            -o -name 'slurm-*.err' \
        \) \
        -delete
}


case "${TARGET}" in
    025deg|1deg)
        validate_run_name "${RUN_NAME}"

        OUTPUT_TARGET="${OUT_ROOT}/${TARGET}/${RUN_NAME}"

        echo "Selected run:"
        echo "  resolution : ${TARGET}"
        echo "  run.name   : ${RUN_NAME}"
        echo "  directory  : ${OUTPUT_TARGET}"

        if [[ ! -d "${OUTPUT_TARGET}" ]]; then
            echo "Run directory does not exist:"
            echo "  ${OUTPUT_TARGET}"
            exit 0
        fi

        echo
        echo "The complete run directory will be removed, including:"
        echo "  NetCDF outputs"
        echo "  run_config.yml"
        echo "  run_manifest.json"
        echo "  run_manifest.rank*.json"
        echo "  validation reports"
        ;;

    all)
        OUTPUT_TARGET="${OUT_ROOT}"

        echo "All model outputs below this directory will be removed:"
        echo "  ${OUTPUT_TARGET}"

        if [[ "${CLEAN_LOGS}" == "1" ]]; then
            echo "SLURM logs will also be removed from:"
            echo "  ${LOG_DIR}"
        fi
        ;;

    *)
        usage
        exit 2
        ;;
esac


echo
read -r -p "Confirm deletion? (y/N): " confirm

if [[ "${confirm}" != "y" && "${confirm}" != "Y" ]]; then
    echo "Cancelled."
    exit 0
fi


if [[ "${TARGET}" == "all" ]]; then
    if [[ -d "${OUT_ROOT}" ]]; then
        find "${OUT_ROOT}" -mindepth 1 -maxdepth 1 -exec rm -rf -- {} +
        echo "Removed all model outputs under:"
        echo "  ${OUT_ROOT}"
    else
        echo "Output directory does not exist:"
        echo "  ${OUT_ROOT}"
    fi
else
    # Safety check: only permit deleting a direct child of
    # Out_ncfile/<resolution>/.
    EXPECTED_PARENT="${OUT_ROOT}/${TARGET}"
    ACTUAL_PARENT="$(cd "$(dirname "${OUTPUT_TARGET}")" && pwd)"

    if [[ "${ACTUAL_PARENT}" != "${EXPECTED_PARENT}" ]]; then
        echo "Safety check failed." >&2
        echo "Expected parent: ${EXPECTED_PARENT}" >&2
        echo "Actual parent:   ${ACTUAL_PARENT}" >&2
        exit 1
    fi

    rm -rf -- "${OUTPUT_TARGET}"

    echo "Removed run directory:"
    echo "  ${OUTPUT_TARGET}"

    # Remove an empty resolution directory after the run is deleted.
    if [[ -d "${EXPECTED_PARENT}" ]] && \
       [[ -z "$(find "${EXPECTED_PARENT}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
        rmdir "${EXPECTED_PARENT}"
    fi
fi


remove_logs

echo "Cleanup complete."