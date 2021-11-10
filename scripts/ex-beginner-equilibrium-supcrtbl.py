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

solution = AqueousPhase(speciate("H O C Na Cl Ca Mg Si"), exclude("organic"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

gases = GaseousPhase("CO2(g) H2O(g)")
gases.setActivityModel(ActivityModelPengRobinson())

minerals = MineralPhases()

system = ChemicalSystem(db, solution, gases, minerals)

state = ChemicalState(system)
state.temperature(60.0, "celsius")
state.pressure(100.0, "bar")
state.set("H2O(aq)"  , 1.0, "kg")
state.set("CO2(g)"   , 1.0, "mol")
state.set("Halite"   , 1.0, "mol")
state.set("Calcite"  , 1.0, "mol")
state.set("Magnesite", 1.0, "mol")
state.set("Quartz"   , 1.0, "mol")

solver = EquilibriumSolver(system)
solver.solve(state)

props = ChemicalProps(state)
props.output("props.txt")

aprops = AqueousProps(state)
aprops.output("aprops.txt")
