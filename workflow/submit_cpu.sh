#!/usr/bin/bash
#SBATCH --job-name=pacemaker-CPU
#SBATCH --account=open-34-9
#SBATCH --partition=qcpu_exp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8 
#SBATCH --cpus-per-task=16

set -e

ml Python/3.9
ml GCCcore/13.2.0

UPLOAD_PATH="$(pwd)/uploads"
LARGE_THRESHOLD=1750
MAX_RUNNING=8
SLEEP_TIME=10

for folder in $(find $UPLOAD_PATH/ -maxdepth 1 -mindepth 1 -type d); do
    FILE="$folder/$(basename "$folder").pckl.gzip"
    FOLDER="$(basename "$folder")"

    if [ ! -f "$FILE" ]; then
        echo "Skipping $FOLDER"
        continue
    fi

    FILE_SIZE=$(stat -c%s "$FILE")
    if [ "$FILE_SIZE" -le "$LARGE_THRESHOLD" ]; then
        cd "$folder"
        while true; do
            if [ "$(pgrep -cx pacemaker)" -ge "$MAX_RUNNING" ]; then
                sleep $SLEEP_TIME
            else
                (OMP_NUM_THREADS=16 pacemaker input.yaml) &
                break
            fi
        done
    fi
done
