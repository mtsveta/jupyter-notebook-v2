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

results_folder = 'results-phrqc2-fugacity-fixed'
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

T = 25.0 # temperature in celsius
P = 1.0  # pressure in bar

solver = EquilibriumSolver(specs)

conditions = EquilibriumConditions(specs)

state = ChemicalState(system)
state.temperature(T, "celsius")
state.set("H2O"        ,  1.0 , "kg")
state.set("Nahcolite"  , 10.00, "mol")
state.set("Natron"     ,  0.00, "mol")
state.set("Trona"      ,  0.00, "mol")
state.set("Na2CO3:H2O" ,  0.00, "mol")
state.set("Na2CO3:7H2O",  0.00, "mol")

def equilibrate(ppCO2, state0):

    conditions.temperature(T, "celsius")
    conditions.pressure(P, "atm")
    conditions.fugacity("CO2", 10**(-ppCO2))

    solver.solve(state, conditions)

    props.update(state)
    aprops.update(state)

    ph = aprops.pH()[0]
    mCO3 = state.speciesAmount("CO3-2")[0]
    mHCO3 = state.speciesAmount("HCO3-")[0]

    return ph, mCO3, mHCO3

num_ppco2s = 71
co2ppressures = np.flip(np.linspace(-5.0, 2.0, num=num_ppco2s))

data_size = 3
data = np.zeros((num_ppco2s, data_size+1))

#print(f'ppCO2   pH      CO3-2   HCO3-')
for i in range(0, num_ppco2s):
    result = equilibrate(co2ppressures[i], state)
    data[i, 0] = co2ppressures[i]
    data[i, 1] = result[0]
    data[i, 2] = result[1]
    data[i, 3] = result[2]
    #print(f'{co2ppressures[i]:4.2f} {result[0]:4.2f} {result[1]:6.4e} {result[2]:6.4e}')

#print(data)
np.savetxt(results_folder + '/m-data.txt', data)
np.savetxt(results_folder + '/m-pH.txt', data[:, 1])
np.savetxt(results_folder + '/m-mCO3.txt', data[:, 2])
np.savetxt(results_folder + '/m-mHCO3.txt', data[:, 3])

import matplotlib.pyplot as plt
colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C7', 'C8', 'C9']

plt.figure()
plt.plot(data[:, 0], data[:, 1], label=f'pH', color=colors[1])
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('pH [-]')
plt.grid()
plt.savefig(results_folder + '/' + 'pH-vs-ppCO2.png', bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(data[:, 0], data[:, 2], label=f'CO3-2', color=colors[2])
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('Amount [mol]')
plt.grid()
plt.savefig(results_folder + '/' + 'mCO32-vs-ppCO2.png', bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(data[:, 0], data[:, 3], label=f'HCO3-', color=colors[3])
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('Amount [mol]')
plt.grid()
plt.savefig(results_folder + '/' + 'mHCO3-vs-ppCO2.png', bbox_inches='tight')
plt.close()

