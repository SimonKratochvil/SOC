#!/usr/bin/bash
#SBATCH --job-name=pace_select-CPU
#SBATCH --account=open-34-9
#SBATCH --partition=qcpu_exp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8 
#SBATCH --cpus-per-task=16


UPLOAD_PATH="$(pwd)/../uploads"
MAX_RUNNING=8
SLEEP_TIME=10

for folder in $(find $UPLOAD_PATH/ -maxdepth 1 -mindepth 1 -type d); do

    FILE="$folder/$(basename "$folder").pckl.gzip"
    POTENTIAL="$folder/output_potential.yaml"
    FOLDER="$(basename "$folder")"

    if [ ! -f "$FILE" ] || [ ! -f "$POTENTIAL" ]; then
        echo "Skipping $FOLDER (missing FILE or POTENTIAL)"
        continue
    fi
    NUM_ENV=$(python3 - <<EOF
import pickle, gzip
with gzip.open("$FILE","rb") as f:
    df = pickle.load(f)
print(sum(len(s) for s in df["ase_atoms"]))
EOF
)

    while [ "$(pgrep -cx pace_select)" -ge "$MAX_RUNNING" ]; do
        sleep "$SLEEP_TIME"
    done
    echo "Launching pace_select in $FOLDER ($NUM_ENV environments)"
    cd "$folder"
    (OMP_NUM_THREADS=16 pace_select -p "$POTENTIAL" "$FILE" &> select.log) &
done

while [ "$(pgrep -cx pace_select)" -gt 0 ]; do
    sleep "$SLEEP_TIME"
done

echo "All jobs done."

