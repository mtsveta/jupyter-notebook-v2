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

db = PhreeqcDatabase("phreeqc.dat")

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

# Define initial equilibrium state
state = ChemicalState(system)
state.setTemperature(T, "celsius")
state.setPressure(P, "bar")
# Take plenty of solution to imitate its equilibration
state.setSpeciesMass("H2O"  , 1.0e-9, "kg")
state.setSpeciesAmount("Na+", 1.2e-9, "mmol")
state.setSpeciesAmount("Cl-", 1.2e-9, "mmol")
# Exchanger site
state.setSpeciesAmount("CaX2", 0.4065, "mol")
state.setSpeciesAmount("KX"  , 0.1871, "mol")

# Define equilibrium solver and equilibrate given initial state with input conditions
solver = EquilibriumSolver(system)
solver.solve(state)
print(state)

aqprops = AqueousProps(state)
print("I  = %f mol/kgw" % aqprops.ionicStrength()[0])
print("pH = %f" % aqprops.pH()[0])
print("pE = %f" % aqprops.pE()[0])

exprops = IonExchangeProps(state)
print(exprops)
