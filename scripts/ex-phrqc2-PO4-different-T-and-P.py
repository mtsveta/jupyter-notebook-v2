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
from math import *

results_folder = 'results-phrqc2-PO4'
os.system('mkdir -p ' + results_folder)

db = PhreeqcDatabase.fromFile('databases/phreeqc-toner-catling.dat') # if running from tutorials folder

print("Database:\n---------------------")
for species in db.species():
    print(species.name())
#input()

# PO4-3
# H2PO4-
# CaH2PO4+
# HPO4-2
# CaHPO4

# FeH2PO4+
# FeH2PO4+2
# FeHCO3+
# FeHPO4
# FeHPO4+
# H3PO4
# KHPO4-
# MgH2PO4+
# MgHCO3+
# MgHPO4
# MgOH+
# MgPO4-

solution = AqueousPhase(speciate("H O C Na Cl Ca P"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

#minerals = MineralPhases("Natron Nahcolite Trona Na2CO3:H2O Na2CO3:7H2O")
minerals = MineralPhases("Na2(HPO4):12H2O Na2(HPO4):7H2O Na2(HPO4):2H2O Na2(HPO4) Na(H2PO4):2H2O Na(H2PO4):H2O Na(H2PO4) "
                         "Natron Nahcolite Trona Na2CO3:H2O Na2CO3:7H2O "
                         "Halite")


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

opts = EquilibriumOptions()
opts.epsilon = 1e-13
solver.setOptions(opts)


def equilibrate_pure(T, ppCO2):

    conditions.temperature(T, "celsius")
    conditions.pressure(1.0, "bar")
    conditions.fugacity("CO2", 10**(ppCO2), 'bar')

    state = ChemicalState(system)
    state.set("H2O", 1.0, "kg")
    state.set("CO2", 100, "mol")
    state.set("Na2(HPO4):7H2O", 10.00, "mol")

    res = solver.solve(state, conditions)

    if not res.optima.succeeded:
        print(f"ERROR: The optimization solver hasn't converged for T = {T} C and ppCO2 = {ppCO2}")
        return nan, nan

    props.update(state)
    aprops.update(state)

    pH = aprops.pH()[0]
    moleP = props.elementAmountInPhase("P", "AqueousPhase")[0]

    #print("pH = ", pH)
    if pH < 5 or pH > 13:
        return nan, nan
    else:
        return pH, moleP

def equilibrate_mixed(T, ppCO2):

    conditions.temperature(T, "celsius")
    conditions.pressure(1.0, "bar")
    conditions.fugacity("CO2", 10**(ppCO2), 'bar')

    state = ChemicalState(system)
    state.set("H2O", 1.0, "kg")
    state.set("CO2", 100, "mol")

    # All the Phosphat species
    state.set("Na2(HPO4):12H2O", 10.00, "mol")
    state.set("Na2(HPO4):7H2O", 10.00, "mol")
    # state.set("Na2(HPO4):2H2O", 10.00, "mol") # solver doesn't converge
    # state.set("Na2(HPO4)", 10.00, "mol")
    # state.set("Na(H2PO4):2H2O", 10.00, "mol")
    # state.set("Na(H2PO4):H2O", 10.00, "mol")
    # state.set("Na(H2PO4)", 10.00, "mol")

    # All the sodium minerals
    state.set("Natron", 10.00, "mol")       # Na2CO3:10H2O
    #state.set("Nahcolite", 10.00, "mol")   # NaHCO3
    #state.set("Trona", 10.00, "mol")       # Na3H(CO3)2:2H2O
    #state.set("Na2CO3:H2O", 10.00, "mol")
    #state.set("Na2CO3:7H2O", 10.00, "mol")
    state.set("Halite", 10.00, "mol")

    res = solver.solve(state, conditions)

    if not res.optima.succeeded:
        print(f"ERROR: The optimization solver hasn't converged for T = {T} C and ppCO2 = {ppCO2}")
        return nan, nan

    props.update(state)
    aprops.update(state)

    pH = aprops.pH()[0]
    moleP = props.elementAmountInPhase("P", "AqueousPhase")[0]

    if pH < 5 or pH > 13:
        return nan, nan
    else:
        return pH, moleP

num_temperatures = 101
num_ppressures = 2
temperatures = np.linspace(0, 50.0, num=num_temperatures)
co2ppressures = np.linspace(-3.5, 0.0, num=num_ppressures)

# data_size = 2
# data0  = np.zeros((num_temperatures, data_size+1))
# data35 = np.zeros((num_temperatures, data_size+1))
#
# for i in range(0, num_temperatures):
#     # ppCO2 = -3.5
#     result = equilibrate_pure(temperatures[i], co2ppressures[0])
#     #print(result)
#     #input()
#     if result != "ERROR":
#         data35[i, 0] = temperatures[i]
#         data35[i, 1] = result[0]
#         data35[i, 2] = result[1]
#
#     # ppCO2 = 0.0
#     result = equilibrate_pure(temperatures[i], co2ppressures[1])
#     #input()
#     if result != "ERROR":
#         data0[i, 0] = temperatures[i]
#         data0[i, 1] = result[0]
#         data0[i, 2] = result[1]
#
# np.savetxt(results_folder + '/data0.txt', data0)
# np.savetxt(results_folder + '/data35.txt', data35)
#
# import matplotlib.pyplot as plt
# colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C7', 'C8', 'C9']
#
# plt.figure()
# plt.plot(temperatures, data0[:, 1], label=f'ppCO2 = 0', color=colors[0])
# plt.plot(temperatures, data35[:, 1], label=f'ppCO2 = -3.5', color=colors[1])
#
# plt.legend(loc="best")
# plt.xlabel('ppCO2')
# plt.ylabel('pH [-]')
# plt.grid()
# plt.savefig(results_folder + '/' + 'pH-vs-ppCO2.png', bbox_inches='tight')
# plt.close()
#
# plt.figure()
# plt.plot(temperatures, data0[:, 2], label=f'ppCO2 = 0', color=colors[2])
# plt.plot(temperatures, data35[:, 2], label=f'ppCO2 = -3.5', color=colors[3])
# plt.yscale('log')
# plt.legend(loc="best")
# plt.xlabel('ppCO2')
# plt.ylabel('Amount of P [mole]')
# plt.grid()
# plt.savefig(results_folder + '/' + 'moleP-vs-ppCO2.png', bbox_inches='tight')
# plt.close()

# ##########################################################

data_size = 2
data0  = np.zeros((num_temperatures, data_size+1))
data35 = np.zeros((num_temperatures, data_size+1))

for i in range(0, num_temperatures):
    # ppCO2 = -3.5
    result = equilibrate_mixed(temperatures[i], co2ppressures[0])
    #print(result)
    #input()
    if result != "ERROR":
        data35[i, 0] = temperatures[i]
        data35[i, 1] = result[0]
        data35[i, 2] = result[1]

    # ppCO2 = 0.0
    result = equilibrate_mixed(temperatures[i], co2ppressures[1])
    #input()
    if result != "ERROR":
        data0[i, 0] = temperatures[i]
        data0[i, 1] = result[0]
        data0[i, 2] = result[1]

np.savetxt(results_folder + '/data0-mixed.txt', data0)
np.savetxt(results_folder + '/data35-mixed.txt', data35)

import matplotlib.pyplot as plt
colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C7', 'C8', 'C9']

plt.figure()
plt.plot(temperatures, data0[:, 1], label=f'ppCO2 = 0', color=colors[0])
plt.plot(temperatures, data35[:, 1], label=f'ppCO2 = -3.5', color=colors[1])

plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('pH [-]')
plt.grid()
plt.savefig(results_folder + '/' + 'pH-vs-ppCO2-mixed.png', bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(temperatures, data0[:, 2], label=f'ppCO2 = 0', color=colors[2])
plt.plot(temperatures, data35[:, 2], label=f'ppCO2 = -3.5', color=colors[3])
plt.yscale('log')
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('Amount of P [mole]')
plt.grid()
plt.savefig(results_folder + '/' + 'moleP-vs-ppCO2-mixed.png', bbox_inches='tight')
plt.close()
