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

results_folder = 'results-stream-model-fixed-fugacity'
os.system('mkdir -p ' + results_folder)

db = PhreeqcDatabase.fromFile('databases/phreeqc-extended.dat') # if running from tutorials folder

solution = AqueousPhase(speciate("H O C Na Cl Ca P"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

minerals = MineralPhases("Fluorapatite Hydroxylapatite Calcite")

system = ChemicalSystem(db, solution, minerals)

props = ChemicalProps(system)
aprops = AqueousProps(system)

specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.fugacity("CO2")

solver = EquilibriumSolver(specs)

conditions = EquilibriumConditions(specs)

state = ChemicalState(system)
state.set("H2O", 1.0, "kg")
state.set("CO2", 100.0, "mmol")
state.set("Calcite", 10.00, "mmol")
state.set("Fluorapatite", 10.00, "mmol")
state.set("Hydroxylapatite", 10.00, "mmol")

def equilibrate(T, ppCO2):

    conditions.temperature(T, "celsius")
    conditions.pressure(1.0, "atm")
    conditions.fugacity("CO2", 10**(ppCO2))

    solver.solve(state, conditions)

    props.update(state)
    aprops.update(state)
    pH = aprops.pH()[0]

    return pH

num_temperatures = 3
num_ppressures = 101
temperatures = np.flip(np.linspace(0, 50.0, num=num_temperatures))
co2ppressures = np.linspace(-4.0, 0.0, num=num_ppressures)

# num_temperatures = 1 #101
# num_ppressures = 1 #106
# temperatures = np.array([50.0])
# co2ppressures = np.array([-4.0])

pHs = np.zeros((num_ppressures, num_temperatures))

# print(temperatures)
# print(co2ppressures)
# input()

p_couter = 0
for ppCO2 in co2ppressures:
    pHs[p_couter, :] = np.array([equilibrate(T, ppCO2) for T in temperatures])
    p_couter += 1

np.savetxt(results_folder + '/pHs.txt', pHs)

import matplotlib.pyplot as plt
colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C7', 'C8', 'C9']

plt.figure()
for i in range(0, num_temperatures):
    plt.plot(co2ppressures, pHs[:, i], label=f'{temperatures[i]} C', color=colors[i])
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('pH [-]')
plt.grid()
plt.savefig(results_folder + '/' + 'pH-vs-ppCO2.png', bbox_inches='tight')
plt.close()


