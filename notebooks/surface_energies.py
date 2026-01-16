from ase.build import bulk, surface, diamond100, diamond111
from ase.visualize import view
from ase.optimize import BFGS
from ase.eos import EquationOfState
from pyace import PyACECalculator
from tqdm.notebook import tqdm
import numpy as np
from quippy.potential import Potential

calc = Potential(init_args='Potential xml_label="GAP_2017_6_17_60_4_3_56_165"',
                                               param_filename='gp_iter6_sparse9k.xml')

dia = bulk('Si', 'diamond', a=5.431, cubic=True)

volumes = []
energies = []

for scale in tqdm(np.linspace(0.98, 1.02, 9)):
    struct = dia.copy()
    struct.set_cell(dia.cell.array*scale, scale_atoms=True)
    struct.calc = calc
    energies.append(struct.get_total_energy())
    volumes.append(struct.get_volume())

eos = EquationOfState(volumes, energies, eos="birchmurnaghan")
V0_fit, E0_fit, B0_fit = eos.fit()
a0_fit = V0_fit**(1/3)
B0_fit *= 160.2

primitive = bulk('Si', 'diamond', a=a0_fit, cubic=True)
supercell = primitive.repeat((3, 3, 15))
supercell.calc = calc
supercell_energy = supercell.get_total_energy()

slab_100 = diamond100('Si', size=(3, 3, 30), a=a0_fit, vacuum=10, orthogonal=True, periodic=True)
slab_111 = diamond111('Si', size=(3, 3, 30), a=a0_fit, vacuum=10, orthogonal=False, periodic=True)

def get_surface_energy(struct):
    struct.calc = calc
    dyn = BFGS(struct)
    dyn.run(fmax=0.001)
    slab_energy = struct.get_total_energy()
    slab_cell = struct.get_cell()
    surface_energy = (slab_energy-supercell_energy/len(supercell)*len(struct))/(2*np.linalg.norm(np.cross(slab_cell[0], slab_cell[1])))*16.02
    return f'Surface energy: {surface_energy}'

print(get_surface_energy(slab_100))
print(get_surface_energy(slab_111))
