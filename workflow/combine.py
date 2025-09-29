import pandas as pd
import glob
import os
import subprocess

files = glob.glob('../uploads/**/selected.pkl.gz', recursive=True)

dfs = []
for f in files:
    df = pd.read_pickle(f)
    dfs.append(df)

combined = pd.concat(dfs, ignore_index=True)

os.makedirs('../combined', exist_ok=True)

data = {
    'energy': combined['energy'].to_numpy(),
    'forces': combined['forces'].to_list(),
    'ase_atoms': combined['ase_atoms'].to_list(),
    'energy_corrected': combined['energy'].to_numpy(),
}

final = pd.DataFrame(data)
final.to_pickle('../combined/combined.pckl.gzip', compression='gzip', protocol=4)

print('Data combined, strating shell script')
subprocess.run(['bash', 'combined_training.sh'])
