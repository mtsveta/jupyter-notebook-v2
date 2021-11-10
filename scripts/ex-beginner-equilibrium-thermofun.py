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

db = ThermoFunDatabase("aq17")

solution = AqueousPhase(speciate("H O C Na Cl"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

gases = GaseousPhase("CO2 H2O")
gases.setActivityModel(ActivityModelPengRobinson())

minerals = MineralPhases()

system = ChemicalSystem(db, solution, gases, minerals)

specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.pH()

solver = EquilibriumSolver(specs)

conditions = EquilibriumConditions(specs)
conditions.temperature(60.0, "celsius")
conditions.pressure(100.0, "bar")
conditions.pH(4.0)
conditions.startWith("H2O@", 1.0, "kg")
conditions.startWith("Na+",  1.0, "mol")
conditions.startWith("Cl-",  1.0, "mol")
conditions.startWith("CO2", 10.0, "mol")

state = ChemicalState(system)

solver.solve(state, conditions)

props = ChemicalProps(state)
props.output("props.txt")

aprops = AqueousProps(state)
aprops.output("aprops.txt")
