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

results_folder = 'results-streammodel-fixed-fugacity-different-ppCO2'
os.system('mkdir -p ' + results_folder)

db = PhreeqcDatabase.fromFile('databases/phreeqc-extended.dat') # if running from tutorials folder

solution = AqueousPhase(speciate("H O C Na Cl Ca P"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))
gases = GaseousPhase("CO2(g)")
gases.setActivityModel(ActivityModelPengRobinson())

mineral = MineralPhase("Calcite Fluorapatite Hydroxylapatite")

system = ChemicalSystem(db, solution, mineral, gases)

aprops = AqueousProps(system)
props = ChemicalProps(system)

specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.fugacity("CO2")

conditions = EquilibriumConditions(specs)

solver = EquilibriumSolver(specs)

opts = EquilibriumOptions()
opts.epsilon = 1e-13
solver.setOptions(opts)

P = 1.0  # pressure in atm

# USE SOLUTION 1;
# EQUILIBRIUM_PHASES 1;
# CO2(g)    -4.00 100;
# Hydroxyapatite 0 10 dissolve_only;
# Fluorapatite 0 10 dissolve_only;
# Calcite 0 10 dissolve_only;
# REACTION_TEMPERATURE 1; 0;
# END

def equilibrate(ppCO2, T):

    conditions.temperature(T, "celsius")
    conditions.pressure(P, "atm")
    conditions.fugacity("CO2", 10 ** (ppCO2), 'atm')

    state = ChemicalState(system)
    state.set("H2O"    , 1.0 , "kg")
    state.set("Calcite", 10.0, "mol")
    state.set("Fluorapatite", 10.0, "mol")
    state.set("Hydroxylapatite", 10.0, "mol")
    state.set("CO2",  100.0, "mol")

    solver.solve(state, conditions)

    aprops.update(state)
    props.update(state)

    pH = aprops.pH()[0]
    molP = props.elementAmountInPhase("P", "AqueousPhase")[0]

    return pH, molP

num_temperatures = 1
num_ppco2s = 101
temperatures = np.array([0, 25, 50])
co2ppressures = np.linspace(-4.0, 0.0, num=num_ppco2s)

data_size = 2
data50 = np.zeros((num_ppco2s, data_size + 1))
data25 = np.zeros((num_ppco2s, data_size + 1))
data0 = np.zeros((num_ppco2s, data_size + 1))

for i in range(0, num_ppco2s):

    result = equilibrate(co2ppressures[i], temperatures[2])
    data50[i, 0] = co2ppressures[i]
    data50[i, 1] = result[0]
    data50[i, 2] = result[1]

    result = equilibrate(co2ppressures[i], temperatures[1])
    data25[i, 0] = co2ppressures[i]
    data25[i, 1] = result[0]
    data25[i, 2] = result[1]

    result = equilibrate(co2ppressures[i], temperatures[0])
    data0[i, 0] = co2ppressures[i]
    data0[i, 1] = result[0]
    data0[i, 2] = result[1]

np.savetxt(results_folder + '/m-data-25.txt', data25)
np.savetxt(results_folder + '/m-pH-25.txt', data25[:, 1])
np.savetxt(results_folder + '/m-mP-25.txt', data25[:, 2])
np.savetxt(results_folder + '/m-data-50.txt', data50)
np.savetxt(results_folder + '/m-pH-50.txt', data50[:, 1])
np.savetxt(results_folder + '/m-mP-50.txt', data50[:, 2])
np.savetxt(results_folder + '/m-data-0.txt', data0)
np.savetxt(results_folder + '/m-pH-0.txt', data0[:, 1])
np.savetxt(results_folder + '/m-mP-0.txt', data0[:, 2])

import matplotlib.pyplot as plt
colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C7', 'C8', 'C9']

plt.figure()
plt.plot(co2ppressures, data50[:, 1], label=f'{temperatures[2]} C', color=colors[2])
plt.plot(co2ppressures, data25[:, 1], label=f'{temperatures[1]} C', color=colors[1])
plt.plot(co2ppressures, data0[:, 1], label=f'{temperatures[0]} C', color=colors[0])
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('pH [-]')
plt.grid()
plt.savefig(results_folder + '/' + 'pH-vs-ppCO2.png', bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(co2ppressures, data50[:, 2], label=f'{temperatures[2]} C', color=colors[2])
plt.plot(co2ppressures, data25[:, 2], label=f'{temperatures[1]} C', color=colors[1])
plt.plot(co2ppressures, data0[:, 2], label=f'{temperatures[0]} C', color=colors[0])
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('Amount of P [mol]')
plt.grid()
plt.savefig(results_folder + '/' + 'mP-vs-ppCO2.png', bbox_inches='tight')
plt.close()


