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

solution = AqueousPhase(speciate("Na Cl C Ca Mg Si"), exclude("organic"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

gases = GaseousPhase("CO2(g) H2O(g)")
gases.setActivityModel(ActivityModelPengRobinson())

minerals = MineralPhases("Halite Calcite Magnesite Dolomite Quartz")

system = ChemicalSystem(db, solution, gases, minerals)

statex = ChemicalState(system)
statex.temperature(60.0, "celsius")
statex.pressure(100.0, "bar")
statex.set("H2O(aq)"  , 1.00, "kg")
statex.set("Halite"   , 1.00, "mol")
statex.set("Calcite"  , 1.00, "mol")
statex.set("Magnesite", 1.00, "mol")
statex.set("Quartz"   , 1.00, "mol")

equilibrate(statex)

propsx = ChemicalProps(statex)
Vx = propsx.volume()
Ux = propsx.internalEnergy()

statex.output("state-expected.txt")
propsx.output("props-expected.txt")

specs = EquilibriumSpecs(system)

idxV = specs.addInput("V")
idxU = specs.addInput("U")

volumeConstraint = ConstraintEquation()
volumeConstraint.id = "VolumeConstraint"
volumeConstraint.fn = lambda props, w: props.volume() - w[idxV]

internalEnergyConstraint = ConstraintEquation()
internalEnergyConstraint.id = "InternalEnergyConstraint"
internalEnergyConstraint.fn = lambda props, w: props.internalEnergy() - w[idxU]

specs.addConstraint(volumeConstraint)
specs.addConstraint(internalEnergyConstraint)

conditions = EquilibriumConditions(specs)
conditions.set("V", Vx)
conditions.set("U", Ux)
conditions.setLowerBoundPressure(1.0, "bar")

state = ChemicalState(system)
state.temperature(25.0, "celsius")
state.pressure(1.0, "bar")
state.set("H2O(aq)"  , 1.00, "kg")
state.set("Halite"   , 1.00, "mol")
state.set("Calcite"  , 1.00, "mol")
state.set("Magnesite", 1.00, "mol")
state.set("Quartz"   , 1.00, "mol")

solver = EquilibriumSolver(specs)

solver.solve(state, conditions)

props = ChemicalProps(state)

state.output("state.txt")
props.output("props.txt")

# Check if props.txt and props-expected.txt are numerically equivalent.
