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

# # Dependence of the pH on added contaminant in water
#
# This tutorial demonstrates how pH is dependent on the added contaminant in the water, affecting the fish life as well
# as the general ecosystem.

# We start by importing the **reaktoro** package:

from reaktoro import *

# To initialize chemical system, we have to start from defining a thermodynamic database and chemical editor (where
# system's phases are defined):

db =  SupcrtDatabase("supcrt98")
aqueousphase = AqueousPhase(speciate("H O Na Cl N"))
aqueousphase.setActivityModel(ActivityModelHKF())

# Construct the chemical system:

system = ChemicalSystem(db, aqueousphase)

for s in system.species():
    print(s.name())

# Below, we initialize chemical state corresponding to the pure water with approximate pH equal to 7:

# +
T = 25.0
P = 1.0

state = ChemicalState(system)
state.setTemperature(T, "celsius")
state.setPressure(P, "bar")
state.setSpeciesMass("H2O(aq)", 1.0, "kg")
# -

# We also define equilibrium solver to be used for range of equilibrium problems:

solver = EquilibriumSolver(system)
solver.solve(state)

# Evaluate pH of the pute water:

aprops = AqueousProps(state)
print("pH of pure water = ", aprops.pH())

# ### Decreasing pH
#
# First, we investigate the behavior of the pH when adding the acidic contaminant to the water.
# We define the auxiliary lists with amounts of acid HCl in the chemical state `hcl_amounts` and corresponding to that
# state list of pH `phs`. Both lists are populated in the loop of 50 steps. We gradually add 0.1 mmol of hydrogen
# chloride and evaluate the pH in the obtained state.

# +
# Initialize lists with HCl amounts and ph values
hcl_initial = 0.0
hcl_delta = 0.1
nsteps = 50
hcl_amounts = [hcl_initial]
phs = [aprops.pH()[0]]

# Run loop of nsteps steps
for i in range(nsteps):

    # Add more hydrogen chlorite to the problem
    state.add("HCl(aq)", hcl_delta, "mmol")

    # Equilibrate state with updated problem
    solver.solve(state)

    # Update aqueous properties
    aprops.update(state)

    # Append new ph
    phs.append(aprops.pH()[0])

    # Append new hydrogen chlorite amount
    hcl_amounts.append(hcl_amounts[-1] + hcl_delta)
# -

# ### Increasing pH
#
# If we add in a chemical contaminant such as ammonia (a compound of nitrogen and hydrogen with the formula NH3,
# colorless gas with a characteristic pungent smell), that can increase the pH and affect fish life.

state = ChemicalState(system)
state.setTemperature(T, "celsius")
state.setPressure(P, "bar")
state.setSpeciesMass("H2O(aq)", 1.0, "kg")
solver.solve(state)
aprops.update(state)

# Define the auxiliary lists with amounts of acid HCl in the chemical state and corresponding to that state pH:

# +
nh3_initial = 0.0
nh3_amounts = [nh3_initial]
phs_increase = [aprops.pH()[0]]

nh3_delta = 0.1
nsteps = 50

for i in range(nsteps):

    # Add more ammonia to the problem
    state.add("NH3(aq)", nh3_delta, "mmol")

    # Equilibrate state with updated problem
    solver.solve(state)

    # Update aqueous properties
    aprops.update(state)

    # Append new ph
    phs_increase.append(aprops.pH()[0])

    # Append new ammonia amount
    nh3_amounts.append(nh3_amounts[-1] + nh3_delta)
# -

# Let us plot pH as a function of the HCl and NH<sup>3</sup>(aq) amounts:

# +
import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(2, 1)

ax1.plot(hcl_amounts, phs, label=f"pH", color='C3')
ax1.legend(loc="best")
ax1.set_title('Dependence of pH on HCl amount')
ax1.grid(True)
ax1.set_ylabel('pH [-]')
ax1.set_xlabel(r'HCl amount [mol]')

ax2.plot(nh3_amounts, phs_increase, label=f"pH", color='C4')
ax2.set_title(r'Dependence of pH on NH$\mathsf{_3}$ amount')
ax2.legend(loc="best")
ax2.grid(True)
ax2.set_ylabel('pH [-]')
ax2.set_xlabel(r'NH$\mathsf{_3}$ amount [mol]')

fig.tight_layout()
fig.savefig('ph-dependence-on-contaminants-in-water.png', bbox_inches='tight')
# -

# As expected, the contaminant like HCl can decrease the pH, affecting the ecosystem, whereas
# ammonia removes H<sup>+</sup> proton from the water to produce ammonium and hydroxide and, therefore,
# increases pH.
