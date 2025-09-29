#!/usr/bin/env python3

import requests
from ase import Atoms
from pint import UnitRegistry

ureg = UnitRegistry()
import pandas as pd
import numpy as np
import argparse
import json

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
    params=dict(username='Kratochvil', password='mT69-nEBr)@$Xm8'),
)
token = response.json()['access_token']

mined_ids = []

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
            {'upload_id': {'all': ['8hRRa1cyTGGVXdB3QGfmRQ']}},
        ],
    },
    'pagination': {'page_size': 1000},
    'required': {
        'include': ['entry_id'],
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

# collect the entry_ids in a loop to get around the pagination limits
while True:
    response = requests.post(
        f'{base_url}/entries/query',
        json=myjson,
        headers={'Authorization': f'Bearer {token}'},
    )
    response_json = response.json()
    # print(json.dumps(response_json, indent=2))

    for data in response_json['data']:
        mined_ids.append(data['entry_id'])

    next_value = response_json['pagination'].get('next_page_after_value')
    if not next_value:
        break
    myjson['pagination']['page_after_value'] = next_value

structures = []
energies = []
forces = []

print('Found {} entries'.format(len(mined_ids)))

# Iterate over the intry ids and try to extract the values we need.
for i, item in enumerate(mined_ids):
    print('Parsing item {}'.format(i))
    first_entry_id = item
    response = requests.post(
        f'{base_url}/entries/{first_entry_id}/archive/query',
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
    try:
        calculations = response_json['data']['archive']['run'][0]['calculation']
    except (KeyError, TypeError):
        print('No calculations detected, skipping to next entry.')
        continue
    systems = response_json['data']['archive']['run'][0]['system']
    if len(calculations) != len(systems):
        print('Number of calculations and systems not consistent.')
        continue

    try:
        rmethod = response_json['data']['archive']['results']['method']
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
            if rmethod['simulation']['precision']['k_line_density'] < args.k * 1e-10:
                print('number of k-points cutoff too low')
                continue
        except:
            pass
    except:
        print('No results method section')
        continue

    for i, c in enumerate(calculations):
        # We don't want a charged system so just skip if we have non-zero charge
        try:
            charge = response_json['data']['archive']['run'][0]['method'][0][
                'electronic'
            ]['charge']
            if charge != 0.0:
                print('Charged system')
                continue
        except (KeyError, TypeError):
            pass

        # forces
        try:
            forces_total = (
                response_json['data']['archive']['run'][0]['calculation'][i]['forces'][
                    'total'
                ]['value']
                * ureg.newton
            )
        except (KeyError, TypeError):
            # print("forces not found")
            try:
                forces_total = (
                    response_json['data']['archive']['run'][0]['calculation'][i][
                        'forces'
                    ]['total']['value_raw']
                    * ureg.newton
                )
            except (KeyError, TypeError):
                try:
                    forces_total = (
                        response_json['data']['archive']['run'][0]['calculation'][i][
                            'forces'
                        ]['free']['value']
                        * ureg.newton
                    )
                except (KeyError, TypeError):
                    print('No forces detected, skipping to next entry.')
                    continue

        # system
        atoms = response_json['data']['archive']['run'][0]['system'][i]['atoms']

        m1 = atoms.get('labels')

        m2 = atoms.get('positions')

        m3 = atoms.get('periodic')

        m4 = atoms.get('lattice_vectors')

        if m1 is None or m2 is None or m3 is None or (m4 is None and m3[0] is True):
            print('Missing data, skipping to the next entry.', m1, m2, m3, m4)
            continue

        labels = response_json['data']['archive']['run'][0]['system'][i]['atoms'][
            'labels'
        ]
        positions = (
            np.array(
                response_json['data']['archive']['run'][0]['system'][i]['atoms'][
                    'positions'
                ]
            )
            * 1e10
        )
        periodic = response_json['data']['archive']['run'][0]['system'][i]['atoms'][
            'periodic'
        ]

        if m4 is not None:
            lattice_vectors = (
                np.array(
                    response_json['data']['archive']['run'][0]['system'][i]['atoms'][
                        'lattice_vectors'
                    ]
                )
                * 1e10
            )
            ase_atoms = Atoms(
                labels, positions=positions, pbc=periodic, cell=lattice_vectors
            )
        else:
            ase_atoms = Atoms(labels, positions=positions, pbc=periodic)

        # energy
        try:
            energy_total = response_json['data']['archive']['run'][0]['calculation'][i][
                'energy'
            ]['total']['value']
        except (KeyError, TypeError):
            print('No energies detected, skipping to the next entry.')
            continue

        if ase_atoms.get_cell().handedness == -1:
            print("Pacemaker can't handle left handed cells, skip it")
            continue

        if len(ase_atoms) != len(forces_total):
            print('Inconsistent number of atoms and forces')
            continue

        energy_total_eV = energy_total * 6.24150907 * (10**18)

        structures.append(ase_atoms)
        energies.append(energy_total_eV)
        forces.append(forces_total.to('eV/angstrom').m)

reference_energy = 0

print('Mined {} configurations'.format(len(energies)))
print('Mined {} atoms'.format(sum([len(i) for i in structures])))

if structures:
    for _ in range(10):
        structures.append(structures[-1])
        energies.append(energies[-1])
        forces.append(forces[-1])

data = {
    'energy': energies,
    'forces': forces,
    'ase_atoms': structures,
    'energy_corrected': energies,
}

df = pd.DataFrame(data)

df.to_pickle('my_test_data.pckl.gzip', compression='gzip', protocol=4)
