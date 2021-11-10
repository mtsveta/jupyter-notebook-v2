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


db = SupcrtDatabase("supcrtbl")

solution = AqueousPhase(speciate("Na Cl C"), exclude("organic"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

gases = GaseousPhase("CO2(g) H2O(g)")
gases.setActivityModel(ActivityModelPengRobinson())

system = ChemicalSystem(db, solution, gases)

specs = EquilibriumSpecs(system)
specs.temperature()
specs.phaseAmount("GaseousPhase")

solver = EquilibriumSolver(specs)

state = ChemicalState(system)
state.setTemperature(50.0, "celsius")
state.setPressure(300.0, "bar")
state.set("H2O(aq)", 1.0, "kg")
state.set("Na+",     1.0, "mol")
state.set("Cl-",     1.0, "mol")
state.set("CO2(aq)", 1.0, "mol")

conditions = EquilibriumConditions(specs)
conditions.temperature(50.0, "celsius")
conditions.phaseAmount("GaseousPhase", 1e-10, "mol")
conditions.setLowerBoundPressure(1.0, "bar")
conditions.setUpperBoundPressure(1000.0, "bar")

solver.solve(state, conditions)

print(f"Pressure: {state.pressure() * 1e-5} bar")
