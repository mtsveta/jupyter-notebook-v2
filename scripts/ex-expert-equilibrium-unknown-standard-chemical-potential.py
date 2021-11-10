# -*- coding: utf-8 -*-
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
from autodiff import abs


T = 60.0 + 273.15 # temperature in K
P = 10.0 * 1e5    # pressure in Pa


def computeSolubilityCaCO3(system):
    state = ChemicalState(system)
    state.setTemperature(T)
    state.setPressure(P)
    state.set("H2O(aq)", 1.0, "kg")
    state.set("Calcite", 1.0, "mol")

    equilibrate(state)

    aprops = AqueousProps(state)

    return aprops.elementMolality("Ca")


db = SupcrtDatabase("supcrtbl")

solution = AqueousPhase(speciate("C Ca"), exclude("organic"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

mineral = MineralPhase("Calcite")

system = ChemicalSystem(db, solution, mineral)

solubilityCaCO3 = computeSolubilityCaCO3(system)

constraint = ConstraintEquation()
constraint.id = "solubility[CaCO3]"

aprops = AqueousProps(system)
def fn(props, w):
    aprops.update(props)
    return aprops.elementMolality("Ca") - solubilityCaCO3

constraint.fn = fn

specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.addUnknownStandardChemicalPotential("Calcite")
specs.addConstraint(constraint)

state = ChemicalState(system)
state.setTemperature(T)
state.setPressure(P)
state.set("H2O(aq)", 1.0, "kg")
state.set("Calcite", 10.0, "mol")

solver = EquilibriumSolver(specs)

solver.solve(state)

G0_calcite_expected = system.species().get("Calcite").props(T, P).G0
G0_calcite_computed = state.equilibrium().p()[0]

print(f"=================================")
print(f"G0(calcite) at 60 °C and 10 bar  ")
print(f"=================================")
print(f"expected: {G0_calcite_expected/1000.0} kJ/mol")
print(f"computed: {G0_calcite_computed/1000.0} kJ/mol")
print(f"   error: {abs((G0_calcite_computed - G0_calcite_expected)/G0_calcite_expected) * 100.0} %")
print(f"=================================")
