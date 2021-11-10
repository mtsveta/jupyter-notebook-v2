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

# Function to calculate equilibrium of the carbonates and seawater
def carbonates_in_seawater(system, solver, T, P):

    # Initial values of the carbonates and mass of water
    n0Calcite = 10.0
    n0Dolomite = 10.0
    water_kg = 1.0

    # Define initial equilibrium state corresponding to the seawater
    state_sw = ChemicalState(system)
    state_sw.setTemperature(T, "celsius")
    state_sw.setPressure(P, "atm")
    state_sw.setSpeciesMass("H2O", 1.0, "kg")
    state_sw.setSpeciesMass("Ca+2", 412.3 * water_kg, "mg")
    state_sw.setSpeciesMass("Mg+2", 1290 * water_kg, "mg")
    state_sw.setSpeciesMass("Na+", 10768.0 * water_kg, "mg")
    state_sw.setSpeciesMass("K+", 399.1 * water_kg, "mg")
    state_sw.setSpeciesMass("Cl-", 19353.0 * water_kg, "mg")
    state_sw.setSpeciesMass("HCO3-", 141.682 * water_kg, "mg")
    state_sw.setSpeciesMass("SO4-2", 2712.0 * water_kg, "mg")

    # Calculate chemical state corresponding to the seawater
    res = solver.solve(state_sw)

    # Throw exception if the equilibrium couldn't be found
    if not res.optima.succeeded:
        raise RuntimeError("Equilibrium calculation did not succeed!")

    # Add carbonates
    state_sw.setSpeciesAmount("Dolomite", n0Dolomite, "mol")
    state_sw.setSpeciesAmount("Calcite", n0Calcite, "mol")

    # Equilibrate the seawater with carbonates
    solver.solve(state_sw)

    # Fetch values of the specified species
    nDolomite = state_sw.speciesAmount("Dolomite")
    nCalcite = state_sw.speciesAmount("Calcite")
    nCa2 = state_sw.speciesAmount("Ca+2")
    nMg2 = state_sw.speciesAmount("Mg+2")
    nH = state_sw.speciesAmount("H+")
    nHCO3 = state_sw.speciesAmount("HCO3-")

    return (nCalcite[0], nDolomite[0], nCa2[0], nMg2[0], nH[0], nHCO3[0])

# Initialize a thermodynamic database
db = PhreeqcDatabase("pitzer.dat")

# Create an aqueous phase automatically selecting all species with provided elements
aqueousphase = AqueousPhase(speciate("H O C Ca Cl Na K Mg S Si"))
aqueousphase.setActivityModel(chain(
    ActivityModelPitzerHMW()
))

# Create carbonates phases
calcitephase = MineralPhase("Calcite")
dolomitephase = MineralPhase("Dolomite")

# Collecting all above-defined phases
phases = Phases(db)
phases.add(aqueousphase)
phases.add(calcitephase)
phases.add(dolomitephase)

# Construct the chemical system
system = ChemicalSystem(phases)

# Create the equilibrium solver
solver = EquilibriumSolver(system)

# Define the range of temperatures and pressure for the equilibrium calculations
T = np.arange(25.0, 91.0, 5.0)
P = 1.0

# Fetch specific species amounts
species_amounts = [carbonates_in_seawater(system, solver, x, P) for x in T]  # [0] is needed to get the value of autodiff.real
mCalcite = [molals[0] for molals in species_amounts]
mDolomite = [molals[1] for molals in species_amounts]
mCa2 = [molals[2] for molals in species_amounts]
mMg2 = [molals[3] for molals in species_amounts]
mH = [molals[4] for molals in species_amounts]
mHCO3 = [molals[5] for molals in species_amounts]

# Output species amount after equilibration for a range of the
print(" --------------------------------------------------------------------------------")
print("  Final species amounts w.r.t. temperatures")
print(" --------------------------------------------------------------------------------")
print("   T     Calcite   Dolomite         Ca++         Mg++           H+        HCO3-")
for i in range(len(T)):
    print(f"{T[i]:4.0f}  {mCalcite[i]:10.4f} {mDolomite[i]:10.4f} {mCa2[i]:12.4e} {mMg2[i]:12.4e} {mH[i]:12.4e} {mHCO3[i]:12.4e}")
