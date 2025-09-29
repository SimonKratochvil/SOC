#!/usr/bin/bash
UPLOAD_PATH="$(pwd)/uploads"
LARGE_THRESHOLD=1750

for folder in "$UPLOAD_PATH"/*; do
    if [ -d "$folder" ]; then
        FILE="$folder/$(basename "$folder").pckl.gzip"
        FOLDER="$(basename "$folder")"
	if [ ! -f "$FILE" ]; then
	    echo "Skipping $FOLDER"
	    continue
        fi
        FILE_SIZE=$(stat -c%s "$FILE")
        if [[ "$FILE_SIZE" -gt "$LARGE_THRESHOLD" ]]; then    
	    cd "$folder"
            sbatch <<EOT
#!/usr/bin/bash
#SBATCH --job-name=pacemaker-GPU
#SBATCH --account=open-34-9
#SBATCH --partition=qgpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=16

#SBATCH --time=02:00:00

ml Python/3.9
ml GCCcore/13.2.0
ml cuDNN/8.4.1.50-CUDA-11.7.0
ml CUDA/11.7.0

OMP_NUM_THREADS=16 pacemaker "input.yaml"
EOT
        fi
    fi
done
