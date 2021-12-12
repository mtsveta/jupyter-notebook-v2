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

results_folder = 'results-stream-model-PO4'
os.system('mkdir -p ' + results_folder)

db = PhreeqcDatabase.fromFile('databases/phreeqc-extended.dat') # if running from tutorials folder

solution = AqueousPhase(speciate("H O C Na Cl Ca P"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

minerals = MineralPhases("Fluorapatite Hydroxylapatite Calcite")

system = ChemicalSystem(db, solution, minerals)

# print("Chemical system content:\n---------------------")
# for species in db.species():
#     print(species.name())
# input()

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
    mFluorapatite = state.speciesAmount("Fluorapatite")[0]
    mCalcite = state.speciesAmount("Calcite")[0]

    deltamHydroxylapatite = 10.0 - state.speciesAmount("Hydroxylapatite")[0]
    mNaHPO4 = state.speciesAmount("NaHPO4-")[0]
    mCaHPO4 = state.speciesAmount("CaHPO4")[0]
    mCaPO4 = state.speciesAmount("CaPO4-")[0]
    molalP = aprops.elementMolality("P")[0]
    moleP = props.elementAmountInPhase("P", "AqueousPhase")[0]
    #print(moleP)
    #input()

    return pH, mFluorapatite, mCalcite, \
           deltamHydroxylapatite, \
           mNaHPO4, mCaHPO4, mCaPO4, \
           molalP

num_temperatures = 101
num_ppressures = 2
temperatures = np.linspace(0, 50.0, num=num_temperatures)
co2ppressures = np.linspace(-3.5, 0.0, num=num_ppressures)

data_size = 8
data0  = np.zeros((num_temperatures, data_size+1))
data35 = np.zeros((num_temperatures, data_size+1))

for i in range(0, num_temperatures):
    # ppCO2 = -3.5
    result = equilibrate(temperatures[i], co2ppressures[0])
    #print(result)
    # input()
    data35[i, 0] = temperatures[i]
    data35[i, 1] = result[0]
    data35[i, 2] = result[1]
    data35[i, 3] = result[2]
    data35[i, 4] = result[3]
    data35[i, 5] = result[4]
    data35[i, 6] = result[5]
    data35[i, 7] = result[6]
    data35[i, 8] = result[7]

    # ppCO2 = 0.0
    result = equilibrate(temperatures[i], co2ppressures[1])
    data0[i, 0] = temperatures[i]
    data0[i, 1] = result[0]
    data0[i, 2] = result[1]
    data0[i, 3] = result[2]
    data0[i, 4] = result[3]
    data0[i, 5] = result[4]
    data0[i, 6] = result[5]
    data0[i, 7] = result[6]
    data0[i, 8] = result[7]


np.savetxt(results_folder + '/data0.txt', data0)
np.savetxt(results_folder + '/data35.txt', data35)

import matplotlib.pyplot as plt
colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C7', 'C8', 'C9']

plt.figure()
plt.plot(temperatures, data0[:, 1], label=f'ppCO2 = 0', color=colors[0])
plt.plot(temperatures, data35[:, 1], label=f'ppCO2 = -3.5', color=colors[1])

plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('pH [-]')
plt.grid()
plt.savefig(results_folder + '/' + 'pH-vs-ppCO2.png', bbox_inches='tight')
plt.close()

# plt.figure()
# plt.plot(temperatures, 1e6 * data0[:, 8], label=f'ppCO2 = 0', color=colors[0])
# plt.legend(loc="best")
# plt.xlabel('ppCO2')
# plt.ylabel('Molalily of P [mmolal]')
# plt.grid()
# plt.savefig(results_folder + '/' + 'molalP-vs-ppCO2-0.png', bbox_inches='tight')
# plt.close()
#
# plt.figure()
# plt.plot(temperatures, 1e8 * data35[:, 8], label=f'ppCO2 = -3.5', color=colors[0])
# plt.legend(loc="best")
# plt.xlabel('ppCO2')
# plt.ylabel('Molalily of P [mmolal]')
# plt.grid()
# plt.savefig(results_folder + '/' + 'molalP-vs-ppCO2-35.png', bbox_inches='tight')
# plt.close()


plt.figure()
plt.plot(temperatures, data0[:, 8], label=f'ppCO2 = 0', color=colors[2])
plt.plot(temperatures, data35[:, 8], label=f'ppCO2 = -3.5', color=colors[3])
plt.yscale('log')
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('Molalily of P [molal]')
plt.grid()
plt.savefig(results_folder + '/' + 'molalP-vs-ppCO2.png', bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(temperatures, data0[:, 6], label=f'ppCO2 = 0', color=colors[2])
plt.plot(temperatures, data35[:, 6], label=f'ppCO2 = -3.5', color=colors[3])
plt.yscale('log')
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('Amount of CaHPO4 [mol]')
plt.grid()
plt.savefig(results_folder + '/' + 'mCaHPO4-vs-ppCO2.png', bbox_inches='tight')
plt.close()


plt.figure()
plt.plot(temperatures, data0[:, 4], label=f'ppCO2 = 0', color=colors[2])
plt.plot(temperatures, data35[:, 4], label=f'ppCO2 = -3.5', color=colors[3])
plt.yscale('log')
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel(r'Solubility of Hydroxylapatite [mol]')
plt.grid()
plt.savefig(results_folder + '/' + 'hydroxylapatite-vs-ppCO2.png', bbox_inches='tight')
plt.close()


