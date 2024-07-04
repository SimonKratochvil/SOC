import requests
import json
from ase import Atoms
from pint import UnitRegistry
ureg = UnitRegistry()
import pandas as pd
import numpy as np
import sys
import argparse

parser = argparse.ArgumentParser(description='Mine selected data from Nomad')
parser.add_argument('-c', help='Program name', type=str)
parser.add_argument('-k', help='K line density', type=np.float64)
parser.add_argument('-b', help='Basis set type', type=str)
parser.add_argument('-sp', help='Spin polarized', type=bool)
parser.add_argument('-sk', help='Smearing kind', type=str)
parser.add_argument('-st', help='Scf threshold energy change', type=np.float64)
parser.add_argument('-pc', help='Planewave cutoff', type=np.float64)
parser.add_argument('-ac', help='Apw cutoff', type=np.float64)
args = parser.parse_args()

print(args.c, args.k, args.b, args.sp, args.sk, args.st, args.pc, args.ac)

base_url = 'http://nomad-lab.eu/prod/v1/api/v1'

# Find first 10 Si calculations that are using PBE functional
myjson={
        'query': {
            'results.material.elements_exclusive': {
                'all': ['Si']
            },
            'results.method.simulation.dft.xc_functional_names': {
                'all': ['GGA_C_PBE', 'GGA_X_PBE']
            }
        },
        'pagination': {
            'page_size': 10000
        },
        'required': {
            'include': ['entry_id'],
        }
    }

if args.c:
    myjson['query']['results.method.simulation.program_name'] = {'all': [args.c]}

if args.k:
    myjson['query']['results.method.simulation.precision.k_line_density'] = {'gte': args.k}

if args.b:
    myjson['query']['results.method.simulation.dft.basis_set_type'] = {'all': [args.b]}

if args.sp:
    myjson['query']['results.method.simulation.dft.spin_polarized'] = {'all': [args.sp]}

if args.sk:
    myjson['query']['results.method.simulation.dft.smearing_kind'] = {'all': [args.sk]}

if args.st:
    myjson['query']['results.method.simulation.dft.scf_threshold_energy_change'] = {'lte': args.st}

if args.pc:
    myjson['query']['results.method.simulation.precision.planewave_cutoff'] = {'gte': args.pc}

if args.ac:
    myjson['query']['results.method.simulation.precision.apw_cutoff'] = {'gte': args.ac}

response = requests.post(f'{base_url}/entries/query', json=myjson)
response_json = response.json()

structures = []
energies = []
forces = []

# Iterate over the intry ids and try to extract the values we need.
for item in response_json['data']:
    first_entry_id = item['entry_id']
    response = requests.post(
        f'{base_url}/entries/{first_entry_id}/archive/query',
        json={
            'required': {
                'run': {
                    'system': {
                        'atoms': '*'
                    },
                    'calculation': {
                        'energy': {
                            'total' : '*'
                        },
                        'forces': {
                            'total' : '*'
                        }
                    }
                }
            }
        })
    response_json = response.json()
    try:
        calculations = response_json['data']['archive']['run'][0]['calculation']
    except (KeyError, TypeError):
        print("No calculations detected, skipping to next entry.")
        continue
    systems = response_json['data']['archive']['run'][0]['system']
    if len(calculations) != len(systems):
        print("Number of calculations and systems not consistent.")
        continue

    for i,c in enumerate(calculations):
        # there is no point in taking every step of longer relax or MD trajectories... they will be highly correclated anyway
        # Just take every 50-th (or first and last of shorter)
        if i % 50 != 0 and i != len(calculations) - 1:
            continue

        #forces
        try:
            forces_total = response_json['data']['archive']['run'][0]['calculation'][i]['forces']['total']['value'] * ureg.newton
        except (KeyError, TypeError):
            print("foces not found")
            try:
                forces_total = response_json['data']['archive']['run'][0]['calculation'][i]['forces']['total']['value_raw'] * ureg.newton
            except (KeyError, TypeError):
                print("No forces detected, skipping to next entry.")
                continue

        #system
        atoms = response_json['data']['archive']['run'][0]['system'][i]['atoms']

        m1 = atoms.get('labels')

        m2 = atoms.get('positions')

        m3 = atoms.get('periodic')

        m4 = atoms.get('lattice_vectors')

        if m1 is None or m2 is None or m3 is None or m4 is None:
            print("Missing data, skipping to the next entry.")
            continue

        labels = response_json['data']['archive']['run'][0]['system'][i]['atoms']['labels']
        positions = np.array(response_json['data']['archive']['run'][0]['system'][i]['atoms']['positions']) * 1e10
        periodic = response_json['data']['archive']['run'][0]['system'][i]['atoms']['periodic']

        lattice_vectors = np.array(response_json['data']['archive']['run'][0]['system'][i]['atoms']['lattice_vectors']) * 1e10
        ase_atoms = Atoms(labels, positions=positions, pbc=periodic, cell=lattice_vectors)

        #energy
        try:
            energy_total = response_json['data']['archive']['run'][0]['calculation'][i]['energy']['total']['value']
        except (KeyError, TypeError):
            print("No energies detected, skipping to the neyt entry.")
            continue

        energy_total_eV = energy_total * 6.24150907 * (10**18)

        structures.append(ase_atoms)
        energies.append(energy_total_eV)
        forces.append(forces_total.to("eV/angstrom").m)
        continue
    continue

reference_energy = 0

data = {'energy': energies,
        'forces': forces,
        'ase_atoms': structures,
        'energy_corrected': energies}

df = pd.DataFrame(data)

df.to_pickle('my_test_data.pckl.gzip', compression='gzip', protocol=4)

    # In general one entry can contain more calculations and system settings which we could extract,
    # right now, just look at the first one...

    # Check that we have all the values that we need, specifically:
    # OK - forces in calculations (value or value_raw)
    # OK - system (we need labels and positions, if periodic is [true,true,true] than we also need the lattice_vectors)
    # OK   if periodic is [false, false, false], than it is just a cluster/molecule and no lattice parameters are needed
    # OK - energy total in calculations
    # OK if we are missing something, skip to the next nomad entry

    # Convert:
    # OK energy from joules to eV
    # OK  forces from newtons to eV/Angstrom
    # OK? system needs to be converted to ase Atoms format, see https://wiki.fysik.dtu.dk/ase/ase/atoms.html
    #
    # OK? Iterate over calculations if multiple present


    # when the conversion are done append the converted energy, forces and ase Atoms to energies, forces and structures lists respectivelly

    # output everything to pacemaker format

    #print(json.dumps(response_json, indent=2))
