#!/bin/bash
# Submit a contiguous range of Monte Carlo samples.
#
# Usage:
#   bash server/submit_mc_batch.sh START_SAMPLE END_SAMPLE [MAX_RUNNING]
#
# Examples:
#   # smoke test: baseline + samples 1-2
#   bash server/submit_mc_batch.sh 0 2 12
#
#   # production in manageable batches
#   bash server/submit_mc_batch.sh 1 25 35
#   bash server/submit_mc_batch.sh 26 50 35

set -euo pipefail

START_SAMPLE="${1:?START_SAMPLE is required}"
END_SAMPLE="${2:?END_SAMPLE is required}"
MAX_RUNNING="${3:-35}"
BANDS_TOTAL="${BANDS_TOTAL:-36}"
PROJECT_DIR="${PROJECT_DIR:-/mnt/beegfs/product/lulc0120}"
CODE_DIR="${PROJECT_DIR}/Code"
RUN_CONFIG="${RUN_CONFIG:-${CODE_DIR}/config/run_mc_1deg.yml}"

if (( END_SAMPLE < START_SAMPLE )); then
    echo "END_SAMPLE must be >= START_SAMPLE" >&2
    exit 2
fi

N_SAMPLES=$(( END_SAMPLE - START_SAMPLE + 1 ))
N_TASKS=$(( N_SAMPLES * BANDS_TOTAL ))
LAST_TASK=$(( N_TASKS - 1 ))

cd "${CODE_DIR}"
echo "Submitting samples ${START_SAMPLE}-${END_SAMPLE}: ${N_TASKS} band tasks, max ${MAX_RUNNING} concurrent"

sbatch \
  --array="0-${LAST_TASK}%${MAX_RUNNING}" \
  --export="ALL,MC_SAMPLE_OFFSET=${START_SAMPLE},BANDS_TOTAL=${BANDS_TOTAL},RUN_CONFIG=${RUN_CONFIG}" \
  server/run_mc_band.sbatch
