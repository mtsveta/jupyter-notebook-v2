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

results_folder = 'results-fixed-pressure-different-ppCO2'
os.system('mkdir -p ' + results_folder)

db = PhreeqcDatabase("phreeqc.dat")

solution = AqueousPhase(speciate("H O C Na Cl Ca"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))
gases = GaseousPhase("CO2(g)")
gases.setActivityModel(ActivityModelPengRobinson())

mineral = MineralPhase("Calcite")

system = ChemicalSystem(db, solution, mineral, gases)

def equilibrate(ppCO2, T):

    state = ChemicalState(system)
    state.set("H2O"    , 1.0 , "kg")
    state.set("Calcite", 10.0, "mol")
    state.set("CO2(g)",  100.0, "mol")

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()

    solver = EquilibriumSolver(specs)

    conditions = EquilibriumConditions(specs)
    conditions.temperature(T, "celsius")
    conditions.pressure(10**(ppCO2), "atm")

    solver.solve(state, conditions)

    aprops = AqueousProps(state)

    pH = aprops.pH()[0]

    return pH

num_temperatures = 1
num_ppco2s = 107
temperatures = np.array([50])
co2ppressures = np.linspace(-4.0, 2.0, num=num_ppco2s)

data_size = 1
data = np.zeros((num_ppco2s, data_size + 1))

for i in range(0, num_ppco2s):
    result = equilibrate(co2ppressures[i], temperatures[0])
    data[i, 0] = co2ppressures[i]
    data[i, 1] = result

np.savetxt(results_folder + '/m-data.txt', data)
np.savetxt(results_folder + '/m-pH.txt', data[:, 1])

import matplotlib.pyplot as plt
colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C7', 'C8', 'C9']

plt.figure()
plt.plot(co2ppressures, data[:, 1], label=f'{temperatures[0]} C', color=colors[0])
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('pH [-]')
plt.grid()
plt.savefig(results_folder + '/' + 'pH-vs-ppCO2.png', bbox_inches='tight')
plt.close()


