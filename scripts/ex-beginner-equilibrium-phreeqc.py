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

import os

print(os.path)

from reaktoro import *

import reaktoro

print(reaktoro.__path__[0])

db = PhreeqcDatabase("phreeqc.dat")

solution = AqueousPhase(speciate("H O C Na Cl"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

gases = GaseousPhase("CO2(g)")
gases.setActivityModel(ActivityModelPengRobinson())

system = ChemicalSystem(db, solution, gases)

T = 25.0 # temperature in celsius
P = 1.0  # pressure in bar

state = ChemicalState(system)
state.temperature(T, "celsius")
state.pressure(P, "bar")
state.setSpeciesMass("H2O"     , 1.0 , "kg")
state.setSpeciesAmount("CO2(g)", 10.0, "mol")
state.setSpeciesAmount("Na+"   , 4.00, "mol")
state.setSpeciesAmount("Cl-"   , 4.00, "mol")

solver = EquilibriumSolver(system)
solver.solve(state)

props = ChemicalProps(state)
props.output("props.txt")

aprops = AqueousProps(state)
aprops.output("aprops.txt")
