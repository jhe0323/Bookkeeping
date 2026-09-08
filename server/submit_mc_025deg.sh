#!/bin/bash
set -euo pipefail

# Usage:
#   bash server/submit_mc_025deg.sh START_SAMPLE END_SAMPLE [START_BAND] [END_BAND]
#
# Re-submission is safe: current completed bands are skipped by fingerprint checks.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATE="${SCRIPT_DIR}/run_mc_025deg.sh"
GENERATED="${SCRIPT_DIR}/generated_mc_025deg"

START_SAMPLE="${1:?START_SAMPLE is required}"
END_SAMPLE="${2:?END_SAMPLE is required}"
START_BAND="${3:-0}"
END_BAND="${4:-119}"
SUBMIT_SLEEP="${SUBMIT_SLEEP:-0.2}"

if (( START_SAMPLE < 0 || END_SAMPLE < START_SAMPLE )); then
    echo "Invalid sample range" >&2
    exit 2
fi

if (( START_BAND < 0 || END_BAND > 119 || START_BAND > END_BAND )); then
    echo "Invalid band range; allowed 0-119" >&2
    exit 2
fi

mkdir -p "${GENERATED}"

for SAMPLE_ID in $(seq "${START_SAMPLE}" "${END_SAMPLE}"); do
    SAMPLE_PAD=$(printf "%06d" "${SAMPLE_ID}")

    for BAND_ID in $(seq "${START_BAND}" "${END_BAND}"); do
        BAND_PAD=$(printf "%03d" "${BAND_ID}")
        JOB_SCRIPT="${GENERATED}/run_mc_025deg_s${SAMPLE_PAD}_b${BAND_PAD}.sh"

        sed             -e "s/__SAMPLE_ID__/${SAMPLE_ID}/g"             -e "s/__SAMPLE_PAD__/${SAMPLE_PAD}/g"             -e "s/__BAND_ID__/${BAND_ID}/g"             -e "s/__BAND_PAD__/${BAND_PAD}/g"             "${TEMPLATE}" > "${JOB_SCRIPT}"

        chmod 750 "${JOB_SCRIPT}"
        echo "Submitting sample=${SAMPLE_ID} band=${BAND_ID}"
        dsub -s "${JOB_SCRIPT}"
        sleep "${SUBMIT_SLEEP}"
    done
done

echo "Submission complete."
