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

db = PhreeqcDatabase("phreeqc-rothmund-kornfeld-convention.dat")

# Define an aqueous phase
solution = AqueousPhase("H2O Na+ Cl- H+ OH- K+ Ca+2 Mg+2")
solution.setActivityModel(ActivityModelHKF())

# Define an ion exchange phase
exchange_species = "NaX KX CaX2"
exchange = IonExchangePhase(exchange_species)
exchange.setActivityModel(ActivityModelIonExchangeGainesThomas())

# Create chemical system
system = ChemicalSystem(db, solution, exchange)

T = 25.0 # temperature in celsius
P = 1.0  # pressure in bar

# Define equilibrium solver and equilibrate given initial state with input conditions
solver = EquilibriumSolver(system)

# Define initial equilibrium state
state = ChemicalState(system)
state.setTemperature(T, "celsius")
state.setPressure(P, "bar")
# Take plenty of solution to imitate its equilibration
state.setSpeciesMass("H2O"   , 1.0 , "kg")
state.setSpeciesAmount("K+"  , 0.1 , "mol")
state.setSpeciesAmount("Ca+2"  , 0.05 , "mol")

# Exchanger site
state.setSpeciesAmount("NaX"  , 0.417, "mol")

solver.solve(state)

exprops = IonExchangeProps(state)
print(exprops)
# [exprops.speciesAmount("KX")[0], exprops.speciesAmount("CaX2")[0]]

