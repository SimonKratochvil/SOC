#!/usr/bin/env python3
from ase import Atoms, io
from pyace import PyACECalculator
from ase.optimize import BFGS
from scipy.optimize import curve_fit
import numpy as np
from ase.constraints import UnitCellFilter

n = int(input("Number of cif files: "))

a = 1

# create ASE atoms
for i in range(n):
    at = io.read("CIFS/Si({}).cif".format(a))
    a = a + 1

    # Create calculator attach it to the Atmos
    at.calc = PyACECalculator('output_potential.yaml')

    # Add unit set filter so that we relax positions and cell simultaneously
    atoms = UnitCellFilter(at)

    # Run relaxation so that the maximum force is bellow 0.001eV/Angstrom
    dyn = BFGS(atoms)
    dyn.run(fmax=0.001)

    #energy = at.get_potential_energy()
    #energies = at.get_potential_energies()

    orig_cell = at.get_cell()
    orig_positions = at.get_positions()
    orig_volume = at.get_volume()

    print(orig_cell)

    f=open("EV.txt", "w")

    strains = np.arange(0.98, 1.02, 0.001)
    E = []

    for i in strains:
        at.set_cell(orig_cell, scale_atoms=True)
        at.set_positions(orig_positions)

        cell = orig_cell * i
        at.set_cell(cell, scale_atoms=True)

        energy2 = at.get_potential_energy()
        E.append(energy2)
        f.write("{} {}\n".format(i, energy2))


        continue
    print(E)

#print(cell)

    # relax the strained structure
# dyn = BFGS(at)
# dyn.run(fmax=0.001)
