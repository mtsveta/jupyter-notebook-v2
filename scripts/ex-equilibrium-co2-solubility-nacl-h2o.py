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
import numpy as np

# Function to calculate solubility of CO2 in the NaCl-brine
def solubility_co2(system, solver, T, P, mNaCl):

    # Initial amount of CO2 gas
    n0CO2g = 10.0

    # Define initial chemical state corresponding to the NaCl-brine of the given concentration
    state = ChemicalState(system)
    state.setTemperature(T, "celsius")
    state.setPressure(P, "bar")
    state.setSpeciesMass("H2O", 1.0, "kg")
    state.setSpeciesAmount("CO2(g)", n0CO2g, "mol")
    state.setSpeciesAmount("Na+", mNaCl, "mol")
    state.setSpeciesAmount("Cl-", mNaCl, "mol")

    # Calculate equilibrium state
    res = solver.solve(state)

    # Throw exception if the equilibrium couldn't be found
    if not res.optima.succeeded:
        raise RuntimeError("Equilibrium calculation did not succeed!")

    # Fetch resulting aqueous properties of the chemical state
    aqprops = AqueousProps(state)

    # Return concentration of the carbon in the aqueous phase
    return aqprops.elementMolality("C")[0]

# Initialize a thermodynamic database
db = PhreeqcDatabase("phreeqc.dat")

# Create an aqueous phase automatically selecting all species with provided elements
aqueousphase = AqueousPhase(speciate("H O C Na Cl"))
aqueousphase.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2"),
))

# Create a gaseous phase
gaseousphase = GaseousPhase("CO2(g)")
gaseousphase.setActivityModel(ActivityModelPengRobinson())

# Collecting all above-defined phases
phases = Phases(db)
phases.add(aqueousphase)
phases.add(gaseousphase)

# Construct the chemical system
system = ChemicalSystem(phases)

# Create the equilibrium solver
solver = EquilibriumSolver(system)

# Define the range of temperatures and pressure for the equilibrium calculations
T = np.arange(25.0, 90.0, 5.0)
P = 100.0

# Calculate CO2 solubilities for the range of the temperatures and different brine concentrations
mCO2_1 = [solubility_co2(system, solver, x, P, mNaCl=1.0) for x in T]
mCO2_2 = [solubility_co2(system, solver, x, P, mNaCl=2.0) for x in T]
mCO2_3 = [solubility_co2(system, solver, x, P, mNaCl=4.0) for x in T]

# Output the results
print(" ----------------------------------------------------------------")
print("  CO2 solubilities w.r.t. temperatures")
print(" ----------------------------------------------------------------")
print("   T    1 mol NaCl-brine   2 mol NaCl-brine   4 mol NaCl-brine")
for i in range(len(T)):
    print(f"{T[i]:4.0f}  {mCO2_1[i]:18.4f} {mCO2_2[i]:18.4f} {mCO2_3[i]:18.4f}")
