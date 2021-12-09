# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: notebooks//ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.1
# ---

from reaktoro import *
import numpy as np
import os
import math

results_folder = 'results-phrqc2-fugacity-fixed-different-ppCO2'
os.system('mkdir -p ' + results_folder)

#db = PhreeqcDatabase.fromFile('../databases/phreeqc-toner-catling.dat') # if running from tutorials folder
db = PhreeqcDatabase.fromFile('databases/phreeqc-toner-catling.dat') # if running from tutorials folder

# print("Database content:\n---------------------")
# for species in db.species():
#     print(species.name())

solution = AqueousPhase(speciate("H O C Na Cl"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

minerals = MineralPhases("Natron Nahcolite Trona Na2CO3:H2O Na2CO3:7H2O")

system = ChemicalSystem(db, solution, minerals)

# print("Chemical system content:\n---------------------")
# for species in db.species():
#     print(species.name())

props = ChemicalProps(system)
aprops = AqueousProps(system)

specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.fugacity("CO2")

P = 1.0  # pressure in bar

conditions = EquilibriumConditions(specs)

# USE SOLUTION 1;
# EQUILIBRIUM_PHASES 1;
# CO2(g)   -5 100;
# Natron    0 0;
# Nahcolite 0 10;
# Trona     0 0;
# Na2CO3:H2O 0 0;
# Na2CO3:7H2O 0 0;
# END

solver = EquilibriumSolver(specs)

opts = EquilibriumOptions()
opts.epsilon = 1e-13
solver.setOptions(opts)

def equilibrate(ppCO2, T):

    conditions.temperature(T, "celsius")
    conditions.pressure(P, "atm")
    conditions.fugacity("CO2", 10 ** (ppCO2), "bar")

    state = ChemicalState(system)
    state.set("H2O", 1.0, "kg")
    state.set("Nahcolite", 10.00, "mol")
    state.set("Natron", 0.00, "mol")
    state.set("Trona", 0.00, "mol")
    state.set("Na2CO3:H2O", 0.00, "mol")
    state.set("Na2CO3:7H2O", 0.00, "mol")
    state.set("CO2",        100.00, "mol")

    res = solver.solve(state, conditions)

    if not res.optima.succeeded:
        print(f"The optimization solver hasn't converged for T = {T} C and ppCO2 = {ppCO2}")
        return math.nan, math.nan, math.nan

    props.update(state)
    aprops.update(state)
    ph = aprops.pH()[0]

    mCO3 = state.speciesAmount("CO3-2")[0]
    mHCO3 = state.speciesAmount("HCO3-")[0]

    return ph, mCO3, mHCO3

num_ppco2s = 71
#num_ppco2s = 8
co2ppressures = np.flip(np.linspace(-5.0, 2.0, num=num_ppco2s))

print(co2ppressures)
input()

data_size = 3
data50 = np.zeros((num_ppco2s, data_size+1))
data25 = np.zeros((num_ppco2s, data_size+1))
data0 = np.zeros((num_ppco2s, data_size+1))

for i in range(0, num_ppco2s):
    result = equilibrate(co2ppressures[i], 0)
    data0[i, 0] = co2ppressures[i]
    data0[i, 1] = result[0]
    data0[i, 2] = result[1]
    data0[i, 3] = result[2]

    result = equilibrate(co2ppressures[i], 25)
    data25[i, 0] = co2ppressures[i]
    data25[i, 1] = result[0]
    data25[i, 2] = result[1]
    data25[i, 3] = result[2]

    result = equilibrate(co2ppressures[i], 50)
    data50[i, 0] = co2ppressures[i]
    data50[i, 1] = result[0]
    data50[i, 2] = result[1]
    data50[i, 3] = result[2]

np.savetxt(results_folder + '/m-data25.txt', data25)
np.savetxt(results_folder + '/m-data0.txt', data0)
np.savetxt(results_folder + '/m-data50.txt', data50)

import matplotlib.pyplot as plt
colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C7', 'C8', 'C9']

plt.figure()
plt.plot(data0[:, 0], data0[:, 1], label=f'0 C', color=colors[1])
plt.plot(data0[:, 0], data25[:, 1], label=f'25 C', color=colors[2])
plt.plot(data0[:, 0], data50[:, 1], label=f'50 C', color=colors[3])
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('pH [-]')
plt.grid()
plt.savefig(results_folder + '/' + 'pH-vs-ppCO2.png', bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(data0[:, 0], data0[:, 2], label=f'0 C', color=colors[1])
plt.plot(data0[:, 0], data25[:, 2], label=f'25 C', color=colors[2])
plt.plot(data0[:, 0], data50[:, 2], label=f'50 C', color=colors[3])
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('Amount of CO3-2 [mol]')
plt.grid()
plt.savefig(results_folder + '/' + 'mCO32-vs-ppCO2.png', bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(data0[:, 0], data0[:, 3], label=f'0 C', color=colors[1])
plt.plot(data0[:, 0], data25[:, 3], label=f'25 C', color=colors[2])
plt.plot(data0[:, 0], data50[:, 3], label=f'50 C', color=colors[3])
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('Amount of HCO3- [mol]')
plt.grid()
plt.savefig(results_folder + '/' + 'mHCO3-vs-ppCO2.png', bbox_inches='tight')
plt.close()

