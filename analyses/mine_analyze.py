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
            {
                'results.material.elements_exclusive': {'all': ['Si']},
            },
            {
                'results.method.simulation.dft.xc_functional_names': [
                    'GGA_C_PBE',
                    'GGA_X_PBE',
                ]
            },
            {'not': {'results.method.simulation.dft.xc_functional_names': ['HF_X']}},
            {'not': {'results.material.structural_type': ['atom']}},
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

# analytical variables
# some variables were not set because those categories where not
# represented in data we worked with so this might need to be changed
# must change -pc back to gte
# uploads = []

n_bulk = 0
n_cluster = 0
n_surface = 0

n_cubic = 0
n_hexagonal = 0
n_monoclinic = 0
n_orthorhombic = 0
n_tetragonal = 0
n_triclinic = 0
n_trigonal = 0

k_line_densities = []
n_vasp = 0
n_gpaw = 0
n_exciting = 0
n_openmx = 0
n_k_vasp = 0
n_k_gpaw = 0
n_k_exciting = 0
n_k_openmx = 0

n_scf = 0

sizes = []

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
                'results': {'material': '*', 'method': '*'},
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
        # there is no point in taking every step of longer relax or MD trajectories... they will be highly correclated anyway
        # Just take every 50-th (or first and last if shorter)
        if i % 50 != 0 and i != len(calculations) - 1:
            continue

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

        # analyze section
        # uploads
        # try:
        # upload_id = response_json['metadata']['upload_id']
        # except (KeyError, TypeError):
        # upload_id = None

        # if upload_id not in uploads and upload_id is not None:
        # uploads.append(upload_id)

        # structures
        try:
            structure = response_json['data']['archive']['results']['material'][
                'structural_type'
            ]
        except (KeyError, TypeError):
            structure = None

        if structure is not None:
            structure = structure.lower()
            if structure == 'bulk':
                n_bulk += 1
            elif structure == 'molecule / cluster':
                n_cluster += 1
            elif structure == 'surface':
                n_surface += 1

        # crystal system
        try:
            crystal_system = response_json['data']['archive']['results']['material'][
                'symmetry'
            ]['crystal_system']
        except (KeyError, TypeError):
            crystal_system = None

        if crystal_system is not None:
            crystal_system = crystal_system.lower()
            if crystal_system == 'cubic':
                n_cubic += 1
            elif crystal_system == 'hexagonal':
                n_hexagonal += 1
            elif crystal_system == 'monoclinic':
                n_monoclinic += 1
            elif crystal_system == 'orthorhombic':
                n_orthorhombic += 1
            elif crystal_system == 'tetragonal':
                n_tetragonal += 1
            elif crystal_system == 'triclinic':
                n_triclinic += 1
            elif crystal_system == 'trigonal':
                n_trigonal += 1

        # k line density
        try:
            k_line_density = response_json['data']['archive']['results']['method'][
                'simulation'
            ]['precision']['k_line_density']
        except (KeyError, TypeError):
            k_line_density = None

        try:
            program_name = response_json['data']['archive']['results']['method'][
                'simulation'
            ]['program_name']
        except (KeyError, TypeError):
            program_name = None

        if program_name is not None:
            program_name = program_name.lower()
            if program_name == 'vasp':
                n_vasp += 1
            elif program_name == 'gpaw':
                n_gpaw += 1
            elif program_name == 'exciting':
                n_exciting += 1
            elif program_name == 'openmx':
                n_openmx += 1

            if program_name == 'vasp' and k_line_density is not None:
                n_k_vasp += 1
            elif program_name == 'gpaw' and k_line_density is not None:
                n_k_gpaw += 1
            elif program_name == 'exciting' and k_line_density is not None:
                n_k_exciting += 1
            elif program_name == 'openmx' and k_line_density is not None:
                n_k_openmx += 1

        if k_line_density is not None:
            k_line_densities.append(k_line_density)

        # SCF threshold energy change
        try:
            scf = response_json['data']['archive']['results']['method']['simulation'][
                'dft'
            ]['scf_threshold_energy_change']
        except (KeyError, TypeError):
            scf = None

        if scf is not None:
            n_scf += 1

        # size
        try:
            size = response_json['data']['archive']['results']['material']['topology'][
                'n_atoms'
            ]
        except (KeyError, TypeError):
            size = len(ase_atoms)

        if size is not None:
            sizes.append(size)

reference_energy = 0

print('Mined {} configurations'.format(len(energies)))
print('Mined {} atoms'.format(sum([len(i) for i in structures])))

print(f'Bulk: {n_bulk}\nMolecule / cluster: {n_cluster}\nSurface: {n_surface}')
print(
    f'Cubic: {n_cubic}\nHexagonal: {n_hexagonal}\nMonoclinic: {n_monoclinic}\nOrthorhombic: {n_orthorhombic}\nTetragonal: {n_tetragonal}\nTriclinic: {n_triclinic}\nTrigonal: {n_trigonal}'
)
# print(f"Distinct uploads: {len(uploads)}")
print(
    f'K line density: {len(k_line_densities)}\nVASP: {n_k_vasp}/{n_vasp}\nGPAW: {n_k_gpaw}/{n_gpaw}\nexciting: {n_k_exciting}/{n_exciting}\nOpenMX: {n_k_openmx}/{n_openmx}'
)
print(f'SCF threshold energy change: {n_scf}')

data = {
    'energy': energies,
    'forces': forces,
    'ase_atoms': structures,
    'energy_corrected': energies,
}

df = pd.DataFrame(data)

df.to_pickle('my_test_data.pckl.gzip', compression='gzip', protocol=4)

forces_seriazible = [f.tolist() for f in forces]

with open('forces_list.json', 'w') as file:
    json.dump(forces_seriazible, file)

with open('sizes_list.json', 'w') as file:
    json.dump(sizes, file)
