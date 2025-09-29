import pandas as pd
import glob
import os
import subprocess
import heapq

files = glob.glob('../uploads/**/selected.pkl.gz', recursive=True)

uploads = []

dfs = []
for f in files:
    df = pd.read_pickle(f)
    dfs.append(df)
    uploads.append((f, len(df['ase_atoms'])))

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

top_10 = heapq.nlargest(10, uploads, key=lambda x: x[1])

print(top_10)
