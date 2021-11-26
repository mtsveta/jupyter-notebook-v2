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

results_folder = 'results-stream-model'
os.system('mkdir -p ' + results_folder)

db = PhreeqcDatabase.fromFile('databases/phreeqc-extended.dat') # if running from tutorials folder

# print("Database content:\n---------------------")
# for species in db.species():
#     print(species.name())

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
#state.set("CO2", 100.0, "mmol")
state.set("Calcite", 10.00, "mol")
state.set("Fluorapatite", 10.00, "mol")
state.set("Hydroxylapatite", 10.00, "mol")

def equilibrate(T, ppCO2):

    conditions.temperature(T, "celsius")
    conditions.pressure(1.0, "atm")
    conditions.fugacity("CO2", 10**(ppCO2))

    solver.solve(state, conditions)

    props.update(state)
    aprops.update(state)
    ph = aprops.pH()[0]

    mPO4 = state.speciesAmount("HPO4-2")[0]
    #mPO4 = aprops.speciesMolality("HPO4-2")[0]

    return ph, mPO4

num_temperatures = 101
num_co2s = 106
temperatures = np.flip(np.linspace(0, 50.0, num=num_temperatures))
co2ppressures = np.linspace(-4.1, 0.1, num=num_co2s)

# num_temperatures = 1 #101
# num_co2s = 1 #106
# temperatures = np.array([50.0])
# co2ppressures = np.array([-4.0])

data_size = 2
data = np.zeros((data_size, num_temperatures, num_co2s))

# print(temperatures)
# print(co2ppressures)

p_couter = 0
for ppCO2 in co2ppressures:
    data[:, :, p_couter] = np.array([equilibrate(T, ppCO2) for T in temperatures]).T
    p_couter += 1
pHs = data[0]
mPO4 = data[1]

np.savetxt(results_folder + '/m-pH.txt', pHs)
np.savetxt(results_folder + '/m-mPO4.txt', mPO4)

import matplotlib.pyplot as plt
colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C7', 'C8', 'C9']

plt.figure()
plt.plot(co2ppressures, pHs[0, :], label=f'50 C', color=colors[1])
plt.plot(co2ppressures, pHs[-51, :], label=f'25 C', color=colors[2])
plt.plot(co2ppressures, pHs[-1, :], label=f'0 C', color=colors[3])
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('pH [-]')
plt.grid()
plt.savefig(results_folder + '/' + 'pH-vs-ppCO2.png', bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(temperatures, mPO4[:, -3], label=f'log10 pCO2 = 0.02', color=colors[1])
plt.plot(temperatures, mPO4[:, 16], label=f'log10 pCO2 = -3.5', color=colors[2])
plt.legend(loc="best")
plt.xlabel('T [degC]')
plt.ylabel('m (HPO4) [molal]')
plt.grid()
plt.savefig(results_folder + '/' + 'mHPO4-vs-T.png', bbox_inches='tight')
plt.close()


