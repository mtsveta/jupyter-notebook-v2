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
import numpy as npy

db_bl = SupcrtDatabase("supcrtbl")
db_98 = SupcrtDatabase("supcrt98")

H2Og_bl = db_bl.species().get("H2O(g)")
H2Og_98 = db_98.species().get("H2O(g)")

CO2g_bl = db_bl.species().get("CO2(g)")
CO2g_98 = db_98.species().get("CO2(g)")

H2Sg_bl = db_bl.species().get("H2S(g)")
H2Sg_98 = db_98.species().get("H2S(g)")

P = 80.0

data = []
for T in npy.arange(25.0, 105.0, 5.0):
    G0_H2Og_bl = H2Og_bl.props(T+273.15, P*1e5).G0[0]/1000
    G0_H2Og_98 = H2Og_98.props(T+273.15, P*1e5).G0[0]/1000
    G0_CO2g_bl = CO2g_bl.props(T+273.15, P*1e5).G0[0]/1000
    G0_CO2g_98 = CO2g_98.props(T+273.15, P*1e5).G0[0]/1000
    G0_H2Sg_bl = H2Sg_bl.props(T+273.15, P*1e5).G0[0]/1000
    G0_H2Sg_98 = H2Sg_98.props(T+273.15, P*1e5).G0[0]/1000
    data.append([T, G0_H2Og_bl, G0_H2Og_98, G0_CO2g_bl, G0_CO2g_98, G0_H2Sg_bl, G0_H2Sg_98])

data = npy.array(data)

print(data)

plt.plot(data[:, 0], data[:, 1], label=f"G0 (kJ/mol) of H2O(g) SUPCRTBL")
plt.plot(data[:, 0], data[:, 2], label=f"G0 (kJ/mol) of H2O(g) SUPCRT98")
plt.legend()
plt.savefig("H2Og.pdf")
plt.close()

plt.plot(data[:, 0], data[:, 3], label=f"G0 (kJ/mol) of CO2(g) SUPCRTBL")
plt.plot(data[:, 0], data[:, 4], label=f"G0 (kJ/mol) of CO2(g) SUPCRT98")
plt.legend()
plt.savefig("CO2g.pdf")
plt.close()

plt.plot(data[:, 0], data[:, 5], label=f"G0 (kJ/mol) of H2S(g) SUPCRTBL")
plt.plot(data[:, 0], data[:, 6], label=f"G0 (kJ/mol) of H2S(g) SUPCRT98")
plt.legend()
plt.savefig("H2Sg.pdf")
plt.close()
