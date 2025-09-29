#!/usr/bin/bash
#SBATCH --job-name=pacemaker-gpu
#SBATCH --account open-34-9
#SBATCH --partition qgpu
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --gpus 1
#SBATCH --time 24:00:00

set -e

ml Python/3.9
ml GCCcore/13.2.0
ml cuDNN/8.4.1.50-CUDA-11.7.0
ml CUDA/11.7.0

python combine_troubleshoot.py
