#!/usr/bin/env python3

from mp_api.client import MPRester
from pymatgen.io.ase import AseAtomsAdaptor
from ase import Atoms, io, Atom
from pyace import PyACECalculator
#from ase.filters import UnitCellFilter
from ase.optimize import BFGS
import numpy as np
import statistics
import math

calc = PyACECalculator('output_potential.yaml')
have_active_set = False
try:
    calc.set_active_set("output_potential.asi")
    have_active_set = True
except:
    pass

# add all other tasks Caclulations -> GGA Structure Optimization -> pick the latest one and add mp code below...
    # we have 40 materials project calculations
with MPRester(api_key="ByjiPidDAqBLRBjr0gbAiYq1SOIL03u5") as mpr:
    tasks_doc = mpr.materials.tasks.search(
        ["mp-2683378",
         "mp-2678018",
         "mp-1261619",
         "mp-1203790",
         "mp-2419299",
         "mp-2419526",
         "mp-1196961",
         "mp-1095376",
         "mp-2683948",
         "mp-1201492",
         "mp-1199894",
         "mp-1204627",
         "mp-1200830",
         "mp-1202745",
         "mp-2620691",
         "mp-1204046",
         "mp-2675050",
         "mp-1249703",
         "mp-2424024",
         "mp-2529013",
         "mp-2360947",
         "mp-1245041",
         "mp-1244971",
         "mp-1245242",
         "mp-1244933",
         "mp-1268191",
         "mp-1244990",
         "mp-2615345",
         "mp-2627750",
         "mp-2459501",
         "mp-1080011",
         "mp-676011",
         "mp-2351841",
         "mp-1120747",
         "mp-1056579",
         "mp-1094057",
         "mp-1120444",
         "mp-2390380",
         "mp-2351833",
         "mp-2350974",
         ]
#       fields=["task_id", "orig_inputs", "calcs_reversed", "output", "last_updated"]
    )
    ace_energies = []
    mp_energies = []
    ids = []
    uncertainties = []

    for task in tasks_doc:
        #print(task.output)
        # Covert MP structure to ASE format
        struct_init = AseAtomsAdaptor.get_atoms(task.output.structure)
        # Attach it to the Atmos
        struct_init.calc = calc
        # Here we could relax the structure, but skip this for now
        # struct = UnitCellFilter(struct_init)
        # Evaluate properties
        # dyn = BFGS(struct)
        # dyn.run(fmax=0.001, steps=1000)
        # If we had an active set, we could check the extrapolation grade here
        # print("final gamma:", max(calc.results['gamma']))

        n_atoms = len(struct_init)

        final_energy_ace = struct_init.get_potential_energy()
        if have_active_set:
            uncertainties.append(max(calc.results['gamma']))

        ace_energies.append(final_energy_ace / n_atoms)
        mp_energies.append(task.output.energy / n_atoms)
        ids.append(task.task_id)

    min_val_ace = min(ace_energies)
    min_val_mp = min(mp_energies)

    ace_energies_final = np.array(ace_energies) - min_val_ace
    mp_energies_final = np.array(mp_energies) - min_val_mp

    if len(ace_energies) != len(mp_energies):
        print("Not consistent number of energies!")

    final_stats = np.array(ace_energies_final) - np.array(mp_energies_final)

    for i,dif in enumerate(final_stats):
        if have_active_set:
            print(f"{ids[i]}: mp above hull: {mp_energies_final[i]}, ace abode hull: {ace_energies_final[i]}, diff: {dif}, gamma: {uncertainties[i]}")
        else:
            print(f"{ids[i]}: mp above hull: {mp_energies_final[i]}, ace abode hull: {ace_energies_final[i]}, diff: {dif}")

    #print(final_stats)

    # average
    num_calc = len(tasks_doc)
    sum = np.sum(final_stats)
    average = sum / num_calc
    print(f"Average: {average}")

    # RMSD
    num_val = len(ace_energies_final)
    power = (np.array(ace_energies_final) - np.array(mp_energies_final))**2
    add = np.sum(power)
    divide = add / num_val
    rmsd = math.sqrt(divide)
    print(f"RMSD: {rmsd}")

    # Average absolute difference
    abs_diff = abs(np.array(ace_energies_final) - np.array(mp_energies_final))
    avg_abs_diff = np.sum(abs_diff) / len(abs_diff)
    print(f"Average absolute difference: \n {avg_abs_diff}")

    # Median of absolute difference
    abs_diff_median = statistics.median(abs_diff)
    print(f"Median of absolute difference: {abs_diff_median}")

    # Max absolute difference
    max_abs_diff = max(abs_diff)
    print(f"Max absolute difference: {max_abs_diff}")

# Tasks
# - normalize energie to get energy per atom
# - subtract the lowest value
# - calculate the sum of squared differences

#    print(mp_energies)
#    mp_energies = np.array(mp_energies) - min(mp_energies)
#    print(mp_energies)
#    print(ace_energies)
#    ace_energies = ace_energies - min(ace_energies)
#    print(ace_energies)
