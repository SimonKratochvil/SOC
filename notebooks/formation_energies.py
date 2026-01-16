from ase.visualize import view
from ase.optimize import BFGS
from ase import Atom
from pyace import PyACECalculator
import numpy as np
from ase.io import read
from quippy.potential import Potential

calc = Potential(init_args='Potential xml_label="GAP_2017_6_17_60_4_3_56_165"',
                                               param_filename='gp_iter6_sparse9k.xml')

perfect = read('Si.vasp')
perfect.calc = calc
perfect_energy = perfect.get_total_energy()

vacancy = perfect.copy()
del vacancy[404]

dumbbell = perfect.copy()
orig_pos = dumbbell[404].position.copy()
d = 2.0
direction = np.array([1.0, 1.0, 0.0])
direction /= np.linalg.norm(direction)
pos1 = orig_pos + (d/2) * direction
pos2 = orig_pos - (d/2) * direction
del dumbbell[404]
dumbbell.append(Atom('Si', position=pos1))
dumbbell.append(Atom('Si', position=pos2))

tetrahedral = read('Si-tetra.vasp')

hexagonal = read('Si-hexa.vasp')

def get_formation_energy(struct, i):
    struct.calc = calc
    dyn = BFGS(struct)
    dyn.run(fmax=0.001)
    formation_energy = struct.get_total_energy()-perfect_energy-i*perfect_energy/len(perfect) 
    return f'Formation energy: {formation_energy}'

print(get_formation_energy(vacancy, -1))
print(get_formation_energy(dumbbell, 1))
print(get_formation_energy(hexagonal, 1))
print(get_formation_energy(tetrahedral, 1))
