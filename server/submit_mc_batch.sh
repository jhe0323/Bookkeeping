#!/bin/bash
set -euo pipefail

# Usage:
#   bash server/submit_mc_batch.sh RESOLUTION START_SAMPLE END_SAMPLE [START_BAND] [END_BAND]

RESOLUTION="${1:?RESOLUTION is required: 1deg or 025deg}"
START_SAMPLE="${2:?START_SAMPLE is required}"
END_SAMPLE="${3:?END_SAMPLE is required}"
START_BAND="${4:-}"
END_BAND="${5:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

case "${RESOLUTION}" in
    1deg)
        TARGET="${SCRIPT_DIR}/submit_mc_1deg.sh"
        ;;
    025deg)
        TARGET="${SCRIPT_DIR}/submit_mc_025deg.sh"
        ;;
    *)
        echo "RESOLUTION must be 1deg or 025deg" >&2
        exit 2
        ;;
esac

if [[ -n "${START_BAND}" && -n "${END_BAND}" ]]; then
    exec bash "${TARGET}" "${START_SAMPLE}" "${END_SAMPLE}" "${START_BAND}" "${END_BAND}"
else
    exec bash "${TARGET}" "${START_SAMPLE}" "${END_SAMPLE}"
fi
