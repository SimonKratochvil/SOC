import numpy as np
from ase import Atoms
from ase.eos import EquationOfState
from ase.visualize import view
from ase.io import read
from ase.build import bulk
from pyace import PyACECalculator
from ase.stress import full_3x3_to_voigt_6_stress, voigt_6_to_full_3x3_stress, full_3x3_to_voigt_6_strain, voigt_6_to_full_3x3_strain
from tqdm.notebook import tqdm
from sympy import Matrix, symbols, simplify, S, linsolve
import pandas as pd
import os

os.environ['OMP_NUM_THREADS'] = '1'

#EOS and BULK MODULUS
diamond_struct = bulk('Si', 'diamond', a=5.431, cubic=True)

volumes = []
energies = []

calc = PyACECalculator('Si_npj_CompMat2021.ace')
for scale in tqdm(np.linspace(0.98, 1.02, 9)):
    struct = diamond_struct.copy()
    struct.set_cell(diamond_struct.cell.array*scale, scale_atoms=True)
    struct.calc = calc
    energies.append(struct.get_total_energy())
    volumes.append(struct.get_volume())

eos = EquationOfState(volumes, energies, eos="birchmurnaghan")
V0_fit, E0_fit, B0_fit = eos.fit()
a0_fit = V0_fit**(1/3)
B0_fit *= 160.2

print('--------------INITIAL----------------')
print(f'lattice parameter is {a0_fit:.4f}Å')
print(f'bulk modulus is {B0_fit:.0f}GPa')
#print(f'minimum energy is {E0_fit/2:.3f}eV/at.')
print('-------------------------------------')

diamond_struct = bulk('Si', 'diamond', a=a0_fit, cubic=True)

#STRAIN-STRESS METHOD
magnitude = 0.02
voigt_strain = magnitude*np.array([ 1, 0, 0, 1, 0, 0])
strain_tensor = voigt_6_to_full_3x3_strain(voigt_strain)

def apply_strain(strain_3x3_tensor, initial):
    new_cell = strain_3x3_tensor @ initial.cell.array
    deformed = initial.copy()
    deformed.set_cell(new_cell, scale_atoms=True)
    return deformed

deformed = apply_strain(strain_tensor, diamond_struct)
deformed.get_volume()/diamond_struct.get_volume()

deformed.calc = calc
stress = deformed.get_stress() * 160.22 # convert to eV/A^3 to GPa
diamond_struct.calc=calc
strain = np.dot(deformed.cell.array, np.linalg.inv(diamond_struct.cell.array))
strain = full_3x3_to_voigt_6_strain(strain)

print('------------STRAIN-STRESS------------')
print(f'B = {(1/3)*(stress[0]/strain[0] + 2*stress[1]/strain[0]):.1f} GPa')
print(f'C11 = {stress[0]/strain[0]:.1f} GPa')
print(f'C12 = {stress[1]/strain[0]:.1f} GPa')
print(f'C44 = {stress[3]/strain[3]:.1f} GPa')
print('-------------------------------------')

#ULICS
def get_ULICS(max_eps: float = 1.5e-2) -> np.ndarray:
    ULICS = max_eps / 6.0 * np.array([
        [1, -2, 3, -4, 5, -6],
        [2, 1, -5, -6, 4, 3],
        [3, 4, -1, 5, 6, -2],
        [4, -3, 6, 1, -2, 5],
        [5, 6, 2, -3, -1, -4],
        [6, -5, -4, 2, -3, 1]
    ], dtype=np.float64)
    return ULICS

magnitude = 0.02
strains = get_ULICS(max_eps=magnitude)
strains = np.concatenate((strains, -strains))

data = []
for i, strain in enumerate(tqdm(strains)):
    deformed = apply_strain(voigt_6_to_full_3x3_strain(strain), diamond_struct)
    deformed.calc = calc
    stress = deformed.get_stress() * 160.22
    data.append((strain, stress))

epsilon, sigma = list(zip(*data)) 
epsilon, sigma = map(np.array, (epsilon, sigma))
Cij = np.matmul(np.linalg.pinv(epsilon), sigma)
Cij = 0.5*(Cij + Cij.T)

def project_cubic(cij: np.ndarray) -> np.ndarray:
    cij = np.array(cij)
    projected_cij = np.zeros((6, 6))
    projected_cij[0, 0] = np.around(np.mean([cij[i, i] for i in np.arange(3)]), 2)
    for i in np.arange(1, 3):
        projected_cij[i, i] = projected_cij[0, 0]
    projected_cij[3, 3] = np.around(np.mean([cij[i, i] for i in np.arange(3, 6)]), 2)
    for i in np.arange(4, 6):
        projected_cij[i, i] = projected_cij[3, 3]
    projected_cij[1, 0] = np.around(np.mean([cij[1, 0], cij[2, 0], cij[2, 1]]), 2)
    projected_cij[2, 0] = projected_cij[1, 0]
    projected_cij[2, 1] = projected_cij[1, 0]
    projected_cij[0, 1] = projected_cij[1, 0]
    projected_cij[0, 2] = projected_cij[1, 0]
    projected_cij[1, 2] = projected_cij[1, 0]
    return projected_cij

projectedCij = project_cubic(Cij)

print('----------------ULICS---------------')
print(f'B = {(1/3)*projectedCij[0,0] + (2/3)*projectedCij[0,1]:.1f} GPa')
print(f'C11 = Cij[0][0] GPa')
print(f'C12 = Cij[0][1] GPa')
print(f'C44 = Cij[3][3]')
print('-------------------------------------')

#ENERGY-STRAIN METHOD
C11, C12, C44, delta = symbols('C_11 C_12 C_44 delta')
C = Matrix([
    [C11, C12, C12, 0, 0, 0],
    [C12, C11, C12, 0, 0, 0],
    [C12, C12, C11, 0, 0, 0],
    [0, 0, 0, C44, 0, 0],
    [0, 0, 0, 0, C44, 0],
    [0, 0, 0, 0, 0, C44]
])
strain_orth = Matrix([
    [delta, -delta, delta**2/(1-delta**2), 0, 0 ,0]
]).T
strain_mono = Matrix([
    [0, 0, delta**2/(4-delta**2), 0, 0, delta]
]).T

strain_energy = S('1/2')*(C*strain_orth).T*strain_orth
simplify(strain_energy)

strain_energy = S('1/2')*(C*strain_mono).T*strain_mono
simplify(strain_energy)

data = []
for i, delta in enumerate(tqdm(np.linspace(-0.010, 0.010, 9))):
    voigt_strain = np.array([delta, -delta,  delta**2/(1-delta**2), 0, 0, 0])
    strain = voigt_6_to_full_3x3_strain(voigt_strain)
    deformed = apply_strain(strain, diamond_struct)
    deformed.calc = calc
    energy = deformed.get_total_energy()
    data.append(dict(delta=delta, energy=energy))

data = pd.DataFrame(data)

quad_fit = np.polyfit(data['delta'], data['energy'], deg=2)

C11C12 = quad_fit[0]/diamond_struct.cell.volume*160.2

linsolve([C11+2*C12-3*B0_fit, C11-C12-C11C12], C11, C12)

data = []
for i, delta in enumerate(tqdm(np.linspace(-0.010, 0.010, 9))):
    voigt_strain = np.array([ 0, 0,  delta**2/(4-delta**2), 0, 0, delta])
    strain = voigt_6_to_full_3x3_strain(voigt_strain)
    deformed = apply_strain(strain, diamond_struct)
    deformed.calc = calc
    energy = deformed.get_total_energy()
    data.append(dict(delta=delta, energy=energy))

data = pd.DataFrame(data)

quad_fit_C44 = np.polyfit(data['delta'], data['energy'], deg=2)
quad_fit_C44

2*quad_fit_C44[0]/diamond_struct.cell.volume*160.2
print('------------ENERGY-STRAIN------------')
print(f'C11 = projectedCij[0][0] GPa')
print(f'C12 = projectedCij[0][1] GPa')
print(f'C44 = projectedCij[3][3]')
print('-------------------------------------')
