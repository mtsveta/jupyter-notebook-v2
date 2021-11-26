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

mineral = MineralPhase("Calcite")

system = ChemicalSystem(db, solution, mineral)

T = 25.0 # temperature in celsius
P = 1.0  # pressure in bar

state = ChemicalState(system)
state.temperature(T, "celsius")
state.pressure(P, "bar")
state.set("H2O"    , 1.0 , "kg")
#state.set("Na+"    , 4.00, "mol")
#state.set("Cl-"    , 4.00, "mol")
#state.set("Calcite", 10.0, "mol")

specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.fugacity("CO2")

solver = EquilibriumSolver(specs)

conditions = EquilibriumConditions(specs)
conditions.temperature(25.0, "celsius")
conditions.pressure(1.0, "atm")
conditions.fugacity("CO2", 10**(-2))

solver.solve(state, conditions)

props = ChemicalProps(state)
props.output("props.txt")

aprops = AqueousProps(state)
aprops.output("aprops.txt")
