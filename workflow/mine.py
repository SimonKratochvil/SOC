#!/usr/bin/env python3

import requests
from ase import Atoms
from pint import UnitRegistry
import pandas as pd
import numpy as np
import argparse
import os
import shutil
import subprocess
from collections import defaultdict

ureg = UnitRegistry()

parser = argparse.ArgumentParser(description='Mine selected data from Nomad')
parser.add_argument('-c', help='Program name', type=str)
parser.add_argument('-k', help='K line density', type=np.float64)
parser.add_argument('-b', help='Basis set type', type=str)
parser.add_argument('-sp', help='Spin polarized', type=bool)
parser.add_argument('-s', help='Structural type', type=str)
parser.add_argument('-sk', help='Smearing kind', type=str)
parser.add_argument('-st', help='Scf threshold energy change', type=np.float64)
parser.add_argument('-pc', help='Planewave cutoff', type=np.float64)
parser.add_argument('-ac', help='Apw cutoff', type=np.float64)
args = parser.parse_args()

print(args.c, args.k, args.b, args.sp, args.sk, args.st, args.pc, args.ac)

base_url = 'http://rosat.physics.muni.cz/nomad-oasis/api/v1'

response = requests.get(
    'http://rosat.physics.muni.cz/nomad-oasis/api/v1/auth/token',
    params=dict(username=os.environ['NOMAD_USER'], password=os.environ['NOMAD_PASS']),
)

token = response.json()['access_token']

entry_ids = []
uploads = defaultdict(list)

myjson = {
    'query': {
        'and': [
            {'results.material.elements_exclusive': {'all': ['Si']}},
            {
                'results.method.simulation.dft.xc_functional_names': [
                    'GGA_C_PBE',
                    'GGA_X_PBE',
                ]
            },
            {'not': {'results.method.simulation.dft.xc_functional_names': ['HF_X']}},
            {'not': {'results.material.structural_type': ['atom']}},
            {'not': {'authors': ['Materials Project']}},
        ],
    },
    'pagination': {'page_size': 1000},
    'required': {
        'include': ['entry_id', 'upload_id'],
    },
    'owner': 'all',
}

if args.c:
    myjson['query']['results.method.simulation.program_name'] = {'all': [args.c]}

if args.b:
    myjson['query']['results.method.simulation.dft.basis_set_type'] = {'all': [args.b]}

if args.sp:
    myjson['query']['results.method.simulation.dft.spin_polarized'] = {'all': [args.sp]}

if args.sk:
    myjson['query']['results.method.simulation.dft.smearing_kind'] = {'all': [args.sk]}

if args.st:
    myjson['query']['results.method.simulation.dft.scf_threshold_energy_change'] = {
        'lte': args.st
    }

if args.s:
    myjson['query']['results.material.structural_type'] = {'all': [args.s]}

if args.ac:
    myjson['query']['results.method.simulation.precision.apw_cutoff'] = {'gte': args.ac}

while True:
    response = requests.post(
        f'{base_url}/entries/query',
        json=myjson,
        headers={'Authorization': f'Bearer {token}'},
    )
    response_json = response.json()

    for data in response_json['data']:
        entry_ids.append(data['entry_id'])
        uploads[data['upload_id']].append(data['entry_id'])

    next_value = response_json['pagination'].get('next_page_after_value')
    if not next_value:
        break
    myjson['pagination']['page_after_value'] = next_value

entry_to_upload = {
    entry_id: upload_id
    for upload_id, entry_list in uploads.items()
    for entry_id in entry_list
}

print(f'Found {len(entry_ids)} entries')
print(f'{len(entry_to_upload)}')
print(f'Found {len(uploads)} uploads')

for i, entry_id in enumerate(entry_ids):
    print(f'Parsing item {i}')
    response = requests.post(
        f'{base_url}/entries/{entry_id}/archive/query',
        headers={'Authorization': f'Bearer {token}'},
        json={
            'required': {
                'run': {
                    'system': {'atoms': '*'},
                    'method': '*',
                    'calculation': {
                        'energy': {'total': '*'},
                        'forces': '*',
                    },
                },
                'results': {'method': '*'},
            }
        },
    )
    response_json = response.json()
    upload_id = entry_to_upload.get(entry_id)
    if upload_id:
        uploads[upload_id].append(response_json)

uploads = {
    upload_id: [entry for entry in entries if isinstance(entry, dict)]
    for upload_id, entries in uploads.items()
}

upload_path = os.path.join(os.getcwd(), 'uploads')
os.makedirs(upload_path, exist_ok=True)
for filename in os.listdir(upload_path):
    file_path = os.path.join(upload_path, filename)
    if os.path.isdir(file_path):
        shutil.rmtree(file_path)

solo_structures = []
solo_energies = []
solo_forces = []
triggered = False

for upload_id, entries in uploads.items():
    structures = []
    energies = []
    forces = []

    for entry_json in entries:
        archive = entry_json['data']['archive']
        run = archive['run'][0]

        try:
            calculation = run['calculation']
        except (KeyError, TypeError):
            print('No calculations detected')
            continue

        system = run['system']
        if len(calculation) != len(system):
            print('Number of calculations and systems not consistent')
            continue

        try:
            rmethod = archive['results']['method']
            try:
                if (
                    rmethod['simulation']['precision']['planewave_cutoff']
                    < args.pc * 1.602177e-19
                ):
                    print('Planewave cutoff too low')
                    continue
            except:
                pass
            try:
                if (
                    rmethod['simulation']['precision']['k_line_density']
                    < args.k * 1e-10
                ):
                    print('Number of k-points cutoff too low')
                    continue
            except:
                pass
        except:
            print('No results method section')
            continue

        for i, c in enumerate(calculation):
            try:
                charge = run['method'][0]['electronic']['charge']
                if charge != 0.0:
                    print('Charged system')
                    continue
            except (KeyError, TypeError):
                pass

            try:
                forces_total = calculation[i]['forces']['total']['value'] * ureg.newton
            except (KeyError, TypeError):
                try:
                    forces_total = (
                        calculation[i]['forces']['total']['value_raw'] * ureg.newton
                    )
                except (KeyError, TypeError):
                    try:
                        forces_total = (
                            calculation[i]['forces']['free']['value'] * ureg.newton
                        )
                    except (KeyError, TypeError):
                        print('No forces detected')
                        continue

            atoms = system[i]['atoms']

            m1 = atoms.get('labels')
            m2 = atoms.get('positions')
            m3 = atoms.get('periodic')
            m4 = atoms.get('lattice_vectors')

            if m1 is None or m2 is None or m3 is None or (m4 is None and m3[0] is True):
                print('Missing data', m1, m2, m3, m4)
                continue

            labels = system[i]['atoms']['labels']
            positions = np.array(system[i]['atoms']['positions']) * 1e10
            periodic = system[i]['atoms']['periodic']

            if m4 is not None:
                lattice_vectors = np.array(system[i]['atoms']['lattice_vectors']) * 1e10
                ase_atoms = Atoms(
                    labels, positions=positions, pbc=periodic, cell=lattice_vectors
                )
            else:
                ase_atoms = Atoms(labels, positions=positions, pbc=periodic)

            try:
                energy_total = calculation[i]['energy']['total']['value']
            except (KeyError, TypeError):
                print('No energies detected')
                continue

            if ase_atoms.get_cell().handedness == -1:
                print("Pacemaker can't handle left handed cells")
                continue

            if len(ase_atoms) != len(forces_total):
                print('Inconsistent number of atoms and forces')
                continue

            energy_total_eV = energy_total * 6.24150907 * (10**18)

            structures.append(ase_atoms)
            energies.append(energy_total_eV)
            forces.append(forces_total.to('eV/angstrom').m)

    reference_energy = 0

    print(f'Upload {upload_id}:')
    print(f'Mined {len(energies)} configurations')
    print(f'Mined {sum([len(i) for i in structures])} atoms')

    if len(energies) == 0:
        print('No data mined')
        continue

    if len(energies) == 1:
        if not triggered:
            os.makedirs(os.path.join(upload_path, 'single'))

        solo_structures.extend(structures)
        solo_energies.extend(energies)
        solo_forces.extend(forces)
        triggered = True

    if len(energies) > 1:
        os.makedirs(os.path.join(upload_path, upload_id))

        data = {
            'energy': energies,
            'forces': forces,
            'ase_atoms': structures,
            'energy_corrected': energies,
        }

        df = pd.DataFrame(data)
        df.to_pickle(
            f'{upload_path}/{upload_id}/{upload_id}.pckl.gzip',
            compression='gzip',
            protocol=4,
        )

if triggered:
    data = {
        'energy': solo_energies,
        'forces': solo_forces,
        'ase_atoms': solo_structures,
        'energy_corrected': solo_energies,
    }

    df = pd.DataFrame(data)
    df.to_pickle(
        f'{upload_path}/single/single.pckl.gzip', compression='gzip', protocol=4
    )

num_uploads = len(
    [f for f in os.listdir(upload_path) if os.path.isdir(os.path.join(upload_path, f))]
)

print(f'Total number of uploads with more than one entry: {num_uploads}')
print('Mining finished, starting shell scripts')

subprocess.run(['bash', 'prepare_uploads.sh'])
