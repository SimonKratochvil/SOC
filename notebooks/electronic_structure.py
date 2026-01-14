from ase import Atoms
from ase.build import bulk
from ase.visualize import view
from ase.eos import EquationOfState
from pyace import PyACECalculator
from ase.io import read
import numpy as np
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
import os

os.environ['OMP_NUM_THREADS'] = '1'

diamond_struct = bulk('Si', 'diamond', a=5.431, cubic=True)

volumes = []
energies = []

calc = PyACECalculator('output_potential.yaml')
for scale in tqdm(np.linspace(0.95, 1.05, 11)):
    struct = diamond_struct.copy()
    struct.set_cell(diamond_struct.cell.array*scale, scale_atoms=True)
    struct.calc = calc
    energies.append(struct.get_total_energy())
    volumes.append(struct.get_volume())

eos = EquationOfState(volumes, energies, eos="birchmurnaghan")
V0_fit, E0_fit, B0_fit = eos.fit()
a0_fit = (V0_fit*4)**(1/3)

print(f'lattice parameter is {a0_fit:.3f}Å')
print(f'bulk modulus is {B0_fit*160.2:.0f}GPa')
print(f'minimum energy is {E0_fit:.3f}eV')
