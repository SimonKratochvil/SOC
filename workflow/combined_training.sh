#!/usr/bin/bash
set -e
DATA_PATH="../combined"
KAPPA=1
L1=0
L2=0

cd $DATA_PATH
FILE="combined.pckl.gzip"
pacemaker -t << EOF
$FILE


400
6

EOF

sed -i "s/\(kappa: \)[0-9.]\+/\1$KAPPA/; \
            s/\(L1_coeffs: \)[^,}]*/\1$L1/; \
            s/\(L2_coeffs: \)[^,}]*/\1$L2/; \
            s/maxiter: 2000/maxiter: 1000/" input.yaml

sbatch <<EOT
#!/usr/bin/bash
#SBATCH --job-name=pacemaker-GPU
#SBATCH --account=open-34-9
#SBATCH --partition=qgpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=16

#SBATCH --time=24:00:00

ml Python/3.9
ml GCCcore/13.2.0
ml cuDNN/8.4.1.50-CUDA-11.7.0
ml CUDA/11.7.0

OMP_NUM_THREADS=16 pacemaker "input.yaml"
EOT

cd ..
