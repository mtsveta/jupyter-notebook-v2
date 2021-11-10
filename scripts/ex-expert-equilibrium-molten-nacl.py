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
from matplotlib import pyplot as plt
from pathlib import Path


N = 300
T = 1100
Pmin = 1
Pmax = 1.7e3

filepath = Path(__file__).parent.parent/"resources/db-molten-nacl.yaml"

db = Database.fromFile(str(filepath))

liquid = LiquidPhase("NaCl(l)")

gases = GaseousPhase(speciate("Na Cl"))

system = ChemicalSystem(db, liquid, gases)

state = ChemicalState(system)
state.temperature(T, "K")
state.set("NaCl(l)", 1.0, "mol")

data = {}
data["P"]              = []
data["NaCl(l)"]        = []
data["Na(g)"]          = []
data["Cl(g)"]          = []
data["Cl2(g)"]         = []
data["NaCl(g)"]        = []
data["Na2Cl2(g,ring)"] = []
data["Na+(g)"]         = []
data["Cl-(g)"]         = []

for i in range(N):
    P = Pmin + i*(Pmax - Pmin)/(N - 1)
    state.pressure(P, "Pa")
    equilibrate(state)
    data["P"].append(P/1e3)
    data["NaCl(l)"].append(state.speciesAmount("NaCl(l)")[0])
    data["Na(g)"].append(state.speciesAmount("Na(g)")[0])
    data["Cl(g)"].append(state.speciesAmount("Cl(g)")[0])
    data["Cl2(g)"].append(state.speciesAmount("Cl2(g)")[0])
    data["NaCl(g)"].append(state.speciesAmount("NaCl(g)")[0])
    data["Na2Cl2(g,ring)"].append(state.speciesAmount("Na2Cl2(g,ring)")[0])
    data["Na+(g)"].append(state.speciesAmount("Na+(g)")[0])
    data["Cl-(g)"].append(state.speciesAmount("Cl-(g)")[0])

plt.subplots(figsize=(10, 5))
plt.plot(data["P"], data["NaCl(l)"], label="NaCl(l)")
plt.plot(data["P"], data["Na(g)"], label="Na(g)")
plt.plot(data["P"], data["Cl(g)"], label="Cl(g)")
plt.plot(data["P"], data["Cl2(g)"], label="Cl2(g)")
plt.plot(data["P"], data["NaCl(g)"], label="NaCl(g)")
plt.plot(data["P"], data["Na2Cl2(g,ring)"], label="Na2Cl2(g,ring)")
plt.plot(data["P"], data["Na+(g)"], label="Na+(g)")
plt.plot(data["P"], data["Cl-(g)"], label="Cl-(g)")
plt.xlabel("P [kPa]")
plt.ylabel("Amount [mol]")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.grid()
plt.tight_layout()
plt.savefig(f"molten-nacl-{T}K-{Pmin}-{Pmax}Pa.pdf")
