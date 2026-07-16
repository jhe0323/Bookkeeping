#!/bin/bash

BASE_DIR="/mnt/beegfs/product/lulc0120"

echo "About to clean the following："
echo "1) $BASE_DIR/slurm_band*"
echo "2) $BASE_DIR/logs/*"
echo "3) $BASE_DIR/outputs/summary_1deg*"

read -p "Confirm？(y/n): " confirm

if [[ $confirm != "y" ]]; then
    echo "Cancel"
    exit 0
fi

echo "Cleaning..."

# clean slurm log
rm -f ${BASE_DIR}/slurm_band*

# clean logs folder
rm -rf ${BASE_DIR}/logs/*

# clean summary files in outputs
rm -f ${BASE_DIR}/outputs/summary_0*

echo "Done ✅"