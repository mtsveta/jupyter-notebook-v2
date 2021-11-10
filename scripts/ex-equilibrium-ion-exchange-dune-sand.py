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
solution = AqueousPhase(speciate("H O C Ca Na Mg Cl"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

# Define an ion exchange phase
exchange_species = "NaX CaX2 MgX2"
exchange = IonExchangePhase(exchange_species)
exchange.setActivityModel(ActivityModelIonExchangeGainesThomas())

# Create chemical system
system = ChemicalSystem(db, solution, exchange)

T = 25.0 # temperature in celsius
P = 1.0  # pressure in atm

# Define initial equilibrium state
state = ChemicalState(system)
state.setTemperature(T, "celsius")
state.setPressure(P, "bar")
# Scale solution recipe to match the values of the PHREEQC examples
state.setSpeciesMass("H2O"   , 1.e6, "kg")
state.setSpeciesAmount("Na+" , 1.10, "kmol")
state.setSpeciesAmount("Mg+2", 0.48, "kmol")
state.setSpeciesAmount("Ca+2", 1.90, "kmol")
# Set the number of exchange assuming that it is completely occupied by sodium
state.setSpeciesAmount("NaX" , 0.06, "mol")

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
