#!/usr/bin/bash
UPLOAD_PATH="uploads"
KAPPA=1
L1=0
L2=0

for folder in $(find $UPLOAD_PATH/ -maxdepth 1 -mindepth 1 -type d); do
    echo "Processing $folder"
    cd "$folder"
    FILE="$(basename "$folder").pckl.gzip"
    pacemaker -t << EOF
$FILE


400
6

EOF

    sed -i "s/\(kappa: \)[0-9.]\+/\1$KAPPA/; \
            s/\(L1_coeffs: \)[^,}]*/\1$L1/; \
            s/\(L2_coeffs: \)[^,}]*/\1$L2/; \
	    s/maxiter: 2000/maxiter: 1000/" input.yaml
    NUM_ENV=$(python3 - <<EOF
import pickle, gzip
with gzip.open("$FILE","rb") as f:
    df = pickle.load(f)
print(sum(len(s) for s in df["ase_atoms"]))
EOF
    )
    if [ "$NUM_ENV" -lt 400 ]; then
	sed -i "s/\(number_of_functions_per_element: \)[0-9]\+/\1$NUM_ENV/" input.yaml
    fi
    cd ../..
done

echo "All folders prepared, submitting jobs..."
