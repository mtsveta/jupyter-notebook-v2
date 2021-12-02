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

results_folder = 'results-phrqc2-pressure-fixed'
os.system('mkdir -p ' + results_folder)

db = PhreeqcDatabase.fromFile('databases/phreeqc-toner-catling.dat') # if running from tutorials folder

# print("Database content:\n---------------------")
# for species in db.species():
#     print(species.name())

solution = AqueousPhase(speciate("H O C Na Cl"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

gases = GaseousPhase("CO2(g)")
gases.setActivityModel(ActivityModelPengRobinson())

minerals = MineralPhases("Natron Nahcolite Trona Na2CO3:H2O Na2CO3:7H2O")

system = ChemicalSystem(db, solution, minerals, gases)
# print("Chemical system content:\n---------------------")
# for species in db.species():
#     print(species.name())

specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()

T = 25.0 # temperature in celsius
ppCO2 = -2

solver = EquilibriumSolver(specs)

conditions = EquilibriumConditions(specs)
conditions.temperature(T, "celsius")
conditions.pressure(10**(ppCO2), "atm")

state = ChemicalState(system)
state.set("H2O", 1.0, "kg")
state.set("Nahcolite", 10.00, "mol")
state.set("Natron", 0.00, "mol")
state.set("Trona", 0.00, "mol")
state.set("Na2CO3:H2O", 0.00, "mol")
state.set("Na2CO3:7H2O", 0.00, "mol")
state.set("CO2(g)", 100, "mol")

solver.solve(state, conditions)

props = ChemicalProps(state)
aprops = AqueousProps(state)

ph = aprops.pH()[0]

print("pH = ", ph)
