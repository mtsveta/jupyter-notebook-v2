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

db = PhreeqcDatabase.fromFile('../databases/phreeqc-extended.dat') # if running from tutorials folder

print("Database content:\n---------------------")
for species in db.species():
    print(species.name())

solution = AqueousPhase(speciate("H O C Na Cl P"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

gases = GaseousPhase("CO2(g)")
gases.setActivityModel(ActivityModelPengRobinson())

minerals = MineralPhases("Fluorapatite Hydroxylapatite Calcite")

system = ChemicalSystem(db, solution, gases, minerals)

print("System content:\n---------------------")
for species in db.species():
    print(species.name())

solver = EquilibriumSolver(system)

props = ChemicalProps(system)
aprops = AqueousProps(system)

def convertPPCO2Bar(ppCO2):
    return 10**ppCO2 # log10(P) = ppCO2

def equilibrate(T, ppCO2):

    P = convertPPCO2Bar(ppCO2)
    state = ChemicalState(system)
    state.temperature(T, "celsius")
    state.pressure(P, "atm")
    state.setSpeciesMass("H2O"     , 1.0 , "kg")
    state.setSpeciesAmount("CO2(g)", 100.0, "mmol")
    state.setSpeciesAmount("Na+"   , 4.00, "mmol")
    state.setSpeciesAmount("Cl-"   , 4.00, "mmol")
    state.setSpeciesAmount("Calcite"        , 10.00, "mmol")
    state.setSpeciesAmount("Fluorapatite"   , 10.00, "mmol")
    state.setSpeciesAmount("Hydroxylapatite", 10.00, "mmol")

    solver.solve(state)

    props.update(state)
    aprops.update(state)

    # PO4 - 3
    # H2PO4 -
    # CaH2PO4 +
    # CaHCO3 +
    # HPO4 - 2

    mPO4 = state.speciesAmount("HPO4-2")[0]

    return mPO4


num_temperatures = 101
num_co2s = 106
temperatures = np.linspace(0, 50.0, num=num_temperatures)
co2ppressures = np.linspace(-4.1, 0.1, num=num_co2s)

mols_PO4 = np.zeros((num_temperatures, num_co2s))

print(temperatures)
print(co2ppressures)

t_couter = 0
for T in temperatures:
    mols_PO4[t_couter, :] = np.array([equilibrate(T, ppCO2) for ppCO2 in co2ppressures])
    t_couter += 1

np.savetxt(results_folder + '/m-PO4.txt', mols_PO4)

