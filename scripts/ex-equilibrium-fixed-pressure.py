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

T = 25.0        # temperature in celsius
P = 1.0         # pressure in bar
ppCO2 = -4.0    # partial pressure of CO2

state = ChemicalState(system)
state.set("H2O"    , 1.0 , "kg")
state.set("Calcite", 10.0, "mol")
state.set("CO2(g)",  100.0, "mmol")

specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()

solver = EquilibriumSolver(specs)

conditions = EquilibriumConditions(specs)
conditions.temperature(50.0, "celsius")
conditions.pressure(10**(ppCO2), "atm")

solver.solve(state, conditions)

props = ChemicalProps(state)
aprops = AqueousProps(state)

