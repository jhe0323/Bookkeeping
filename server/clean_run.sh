#!/bin/bash
set -euo pipefail

# Clean deterministic LULCC outputs:
#   Out_ncfile/<resolution>/<run.name>/
#
# Logs are retained by default.
# Set CLEAN_LOGS=1 only when log removal is intentional.

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
    echo "Optional:"
    echo "  CLEAN_LOGS=1 bash server/clean_run.sh ..."
}


validate_run_name() {
    local name="$1"
    if [[ -z "${name}" ]]; then
        echo "Error: run.name is required." >&2
        exit 2
    fi
    if [[ "${name}" == "." || "${name}" == ".." || "${name}" == *"/"* ]]; then
        echo "Error: invalid run.name: ${name}" >&2
        exit 2
    fi
}


remove_logs() {
    # FIX: old code deleted logs when CLEAN_LOGS != 1 and then deleted matching
    # logs unconditionally. Now this function returns unless deletion is explicit.
    if [[ "${CLEAN_LOGS}" != "1" ]]; then
        echo "Logs retained (set CLEAN_LOGS=1 to remove them)."
        return
    fi

    if [[ ! -d "${LOG_DIR}" ]]; then
        echo "Log directory does not exist: ${LOG_DIR}"
        return
    fi

    find "${LOG_DIR}" -maxdepth 1 -type f \
        \( \
            -name 'lulcc_*.out' \
            -o -name 'lulcc_*.err' \
            -o -name 'mc_*.out' \
            -o -name 'mc_*.err' \
            -o -name 'env_test_*.out' \
            -o -name 'env_test_*.err' \
        \) \
        -print \
        -delete
}


case "${TARGET}" in
    025deg|1deg)
        validate_run_name "${RUN_NAME}"
        OUTPUT_TARGET="${OUT_ROOT}/${TARGET}/${RUN_NAME}"
        EXPECTED_PARENT="${OUT_ROOT}/${TARGET}"

        echo "Selected run:"
        echo "  resolution : ${TARGET}"
        echo "  run.name   : ${RUN_NAME}"
        echo "  directory  : ${OUTPUT_TARGET}"

        if [[ ! -d "${OUTPUT_TARGET}" ]]; then
            echo "Run directory does not exist: ${OUTPUT_TARGET}"
            remove_logs
            exit 0
        fi
        ;;
    all)
        OUTPUT_TARGET="${OUT_ROOT}"
        echo "All model outputs below this directory will be removed:"
        echo "  ${OUTPUT_TARGET}"
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
        find "${OUT_ROOT}" \
            -mindepth 1 \
            -maxdepth 1 \
            -exec rm -rf -- {} +
    fi
else
    EXPECTED_PARENT_REAL="$(cd "${EXPECTED_PARENT}" && pwd)"
    ACTUAL_PARENT_REAL="$(cd "$(dirname "${OUTPUT_TARGET}")" && pwd)"

    if [[ "${ACTUAL_PARENT_REAL}" != "${EXPECTED_PARENT_REAL}" ]]; then
        echo "Safety check failed." >&2
        exit 1
    fi

    rm -rf -- "${OUTPUT_TARGET}"

    if [[ -d "${EXPECTED_PARENT}" ]] && \
       [[ -z "$(find "${EXPECTED_PARENT}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
        rmdir "${EXPECTED_PARENT}"
    fi
fi

remove_logs
echo "Cleanup complete."
