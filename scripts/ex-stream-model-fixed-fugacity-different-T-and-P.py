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

# USE SOLUTION 1;
# EQUILIBRIUM_PHASES 1;
# CO2(g)    -4 100;
# Hydroxyapatite 0 10 dissolve_only;
# Fluorapatite 0 10 dissolve_only;
# Calcite 0 10 dissolve_only;
# REACTION_TEMPERATURE 1; 50;
# END

state = ChemicalState(system)
state.set("H2O", 1.0, "kg")
state.set("CO2", 100.0, "mol")
state.set("Calcite", 10.00, "mol")
state.set("Fluorapatite", 10.00, "mol")
state.set("Hydroxylapatite", 10.00, "mol")

def equilibrate(T, ppCO2):

    conditions.temperature(T, "celsius")
    conditions.pressure(1.0, "atm")
    conditions.fugacity("CO2", 10**(ppCO2), 'atm')

    solver.solve(state, conditions)

    props.update(state)
    aprops.update(state)
    pH = aprops.pH()[0]

    return pH

num_temperatures = 3
num_ppressures = 101
temperatures = np.linspace(0, 50.0, num=num_temperatures)
co2ppressures = np.linspace(-4.0, 0.0, num=num_ppressures)

pHs = np.zeros((num_ppressures, num_temperatures))

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


