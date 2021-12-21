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
import os

results_folder = 'results-streammodel-fixed-pressure-different-ppCO2'
os.system('mkdir -p ' + results_folder)

db = PhreeqcDatabase.fromFile('databases/phreeqc-extended.dat') # if running from tutorials folder

solution = AqueousPhase(speciate("H O C Na Cl Ca P"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))
gases = GaseousPhase("CO2(g)")
gases.setActivityModel(ActivityModelPengRobinson())

mineral = MineralPhase("Calcite Fluorapatite Hydroxylapatite")

system = ChemicalSystem(db, solution, mineral, gases)

aprops = AqueousProps(system)
props = ChemicalProps(system)

def equilibrate(ppCO2, T):

    state = ChemicalState(system)
    state.set("H2O"    , 1.0 , "kg")
    state.set("Calcite", 10.0, "mol")
    state.set("Fluorapatite", 10.0, "mol")
    state.set("Hydroxylapatite", 10.0, "mol")
    state.set("CO2(g)",  100.0, "mol")

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()

    solver = EquilibriumSolver(specs)

    opts = EquilibriumOptions()
    opts.epsilon = 1e-13
    solver.setOptions(opts)

    conditions = EquilibriumConditions(specs)
    conditions.temperature(T, "celsius")
    conditions.pressure(10**(ppCO2), "atm")

    solver.solve(state, conditions)

    aprops.update(state)
    props.update(state)

    pH = aprops.pH()[0]
    #molP = props.elementAmountInPhase("P", "AqueousPhase")[0]
    molP = aprops.elementMolality("P")[0]

    return pH, molP

num_temperatures = 3
num_ppco2s = 7
temperatures = np.array([0, 25, 50])
co2ppressures = np.linspace(-4.0, 2.0, num=num_ppco2s)

data_size = 2
data50 = np.zeros((num_ppco2s, data_size + 1))
data25 = np.zeros((num_ppco2s, data_size + 1))
data0 = np.zeros((num_ppco2s, data_size + 1))

for i in range(0, num_ppco2s):

    result = equilibrate(co2ppressures[i], temperatures[2])
    data50[i, 0] = co2ppressures[i]
    data50[i, 1] = result[0]
    data50[i, 2] = result[1]

    result = equilibrate(co2ppressures[i], temperatures[1])
    data25[i, 0] = co2ppressures[i]
    data25[i, 1] = result[0]
    data25[i, 2] = result[1]

    result = equilibrate(co2ppressures[i], temperatures[0])
    data0[i, 0] = co2ppressures[i]
    data0[i, 1] = result[0]
    data0[i, 2] = result[1]

np.savetxt(results_folder + '/m-data-25.txt', data25)
np.savetxt(results_folder + '/m-pH-25.txt', data25[:, 1])
np.savetxt(results_folder + '/m-mP-25.txt', data25[:, 2])
np.savetxt(results_folder + '/m-data-50.txt', data50)
np.savetxt(results_folder + '/m-pH-50.txt', data50[:, 1])
np.savetxt(results_folder + '/m-mP-50.txt', data50[:, 2])
np.savetxt(results_folder + '/m-data-0.txt', data0)
np.savetxt(results_folder + '/m-pH-0.txt', data0[:, 1])
np.savetxt(results_folder + '/m-mP-0.txt', data0[:, 2])

import matplotlib.pyplot as plt
colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C7', 'C8', 'C9']

plt.figure()
plt.plot(co2ppressures, data50[:, 1], label=f'{temperatures[2]} C', color=colors[2])
plt.plot(co2ppressures, data25[:, 1], label=f'{temperatures[1]} C', color=colors[1])
plt.plot(co2ppressures, data0[:, 1], label=f'{temperatures[0]} C', color=colors[0])
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('pH [-]')
plt.grid()
plt.savefig(results_folder + '/' + 'pH-vs-ppCO2.png', bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(co2ppressures, data50[:, 2], label=f'{temperatures[2]} C', color=colors[2])
plt.plot(co2ppressures, data25[:, 2], label=f'{temperatures[1]} C', color=colors[1])
plt.plot(co2ppressures, data0[:, 2], label=f'{temperatures[0]} C', color=colors[0])
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('Amount of P [mol]')
plt.grid()
plt.savefig(results_folder + '/' + 'mP-vs-ppCO2.png', bbox_inches='tight')
plt.close()


# ##
# # Comparison to PHREEQC output
# ##
#
data0PHREEQC  = np.loadtxt(results_folder + '/StreamModel_T0_wider_range.sel' , skiprows=2)
data25PHREEQC = np.loadtxt(results_folder + '/StreamModel_T25.sel', skiprows=2)
data50PHREEQC = np.loadtxt(results_folder + '/StreamModel_T50.sel', skiprows=2)

# comparison plots
def line_filled_marker(color):
    return {'color': color, 'markersize': 6, 'markeredgewidth': 1.5 }

plt.figure()
plt.plot(data0[:, 0], data0[:, 1], label=f'Reaktoro', color=colors[2])
plt.plot(data0[:, 0], data25[:, 1], label=f'Reaktoro', color=colors[4])
plt.plot(data0[:, 0], data50[:, 1], label=f'Reaktoro', color=colors[6])
plt.plot(data0PHREEQC[:, 3], data0PHREEQC[:, 0], 'D', **line_filled_marker(colors[2]))
plt.plot([], [], 'D', label='T = 0 C, PHREEQC', **line_filled_marker(colors[2]))
plt.plot(data25PHREEQC[:, 3], data25PHREEQC[:, 0], 'D', **line_filled_marker(colors[4]))
plt.plot([], [], 'D', label='T = 25 C, PHREEQC', **line_filled_marker(colors[4]))
plt.plot(data50PHREEQC[:, 3], data50PHREEQC[:, 0], 'D', **line_filled_marker(colors[6]))
plt.plot([], [], 'D', label='T = 50 C, PHREEQC', **line_filled_marker(colors[6]))
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel(f'pH [-]')
plt.title(f'pH for different temperatures')
plt.grid()
plt.savefig(results_folder + '/' + 'reaktoro-phreeqc-pH-vs-ppC02.png', bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(data0[:, 0], data0[:, 2], label=f'Reaktoro', color=colors[3])
plt.plot(data0[:, 0], data25[:, 2], label=f'Reaktoro', color=colors[5])
plt.plot(data0[:, 0], data50[:, 2], label=f'Reaktoro', color=colors[7])
plt.plot(data0PHREEQC[:, 3], data0PHREEQC[:, 2],'D', **line_filled_marker(colors[3]))
plt.plot([], [], 'D', label='T = 0 C, PHREEQC', **line_filled_marker(colors[3]))
plt.plot(data25PHREEQC[:, 3], data25PHREEQC[:, 2], 'D', **line_filled_marker(colors[5]))
plt.plot([], [], 'D', label='T = 25 C, PHREEQC', **line_filled_marker(colors[5]))
plt.plot(data50PHREEQC[:, 3], data50PHREEQC[:, 2], 'D', **line_filled_marker(colors[7]))
plt.plot([], [], 'D', label='T = 50 C, PHREEQC', **line_filled_marker(colors[7]))
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel(f'Amount of P [mol]')
plt.title(f'Amount of P for different temperatures')
plt.grid()
plt.savefig(results_folder + '/' + 'reaktoro-phreeqc-P-vs-ppC02.png', bbox_inches='tight')
plt.close()

# #########################################################################################
#
# num_temperatures = 51
# num_ppco2s = 101
# temperatures =  np.linspace(0.0, 50.0, num=num_temperatures)
# co2ppressures = np.linspace(-4.0, 0.0, num=num_ppco2s)
#
# data_size = 2
# data_pH = np.zeros((num_temperatures, num_ppco2s))
# data_P = np.zeros((num_temperatures, num_ppco2s))
#
# for i in range(0, num_temperatures):
#     for j in range(0, num_ppco2s):
#
#         result = equilibrate(co2ppressures[j], temperatures[i])
#         data_pH[i, j] = result[0]
#         data_P[i, j] = result[1]
#
# import matplotlib as ml
# import matplotlib.pyplot as plt
# from matplotlib.colors import BoundaryNorm
# from matplotlib.ticker import MaxNLocator
# import numpy as np
#
# levels = MaxNLocator(nbins=15).tick_values(data_pH.min(), data_pH.max())
# cmap = ml.cm.get_cmap('rainbow')
# norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
#
# fig, ax = plt.subplots(1, 1)
# im = plt.pcolormesh(co2ppressures, temperatures, data_pH, cmap=cmap, norm=norm, shading='auto')
# fig.colorbar(im, ax=ax)
# ax.set_title('pH [-]')
# plt.savefig(results_folder + '/' + 'pH-vs-T-ppCO2.png', bbox_inches='tight')
# plt.close()
#
# levels = MaxNLocator(nbins=15).tick_values(data_P.min(), data_P.max())
# cmap = ml.cm.get_cmap('rainbow')
# norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
#
# fig, ax = plt.subplots(1, 1)
# im = plt.pcolormesh(co2ppressures, temperatures, data_P, cmap=cmap, norm=norm, shading='auto')
# fig.colorbar(im, ax=ax)
# ax.set_title('Amount of P [mol]')
# plt.savefig(results_folder + '/' + 'molP-vs-T-ppCO2.png', bbox_inches='tight')
# plt.close()
#
# levels = MaxNLocator(nbins=15).tick_values(np.log10(data_P).min(), np.log10(data_P).max())
# cmap = ml.cm.get_cmap('rainbow')
# norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
#
# fig, ax = plt.subplots(1, 1)
# im = plt.pcolormesh(co2ppressures, temperatures, np.log10(data_P), cmap=cmap, norm=norm,  shading='auto')
# fig.colorbar(im, ax=ax)
# ax.set_title('Amount log10(P) [mol]')
# plt.savefig(results_folder + '/' + 'log10P-vs-T-ppCO2.png', bbox_inches='tight')
# plt.close()
#
# # Plotting with `contourf` function
#
# maps = ["Blues", "Greys",  "PuBu", "Oranges", "YlOrBr", "Reds", "Purples", "PuRd",  "Greens"]
# dx = co2ppressures[1] - co2ppressures[0]
# dy = temperatures[1] - temperatures[0]
# y, x = np.mgrid[slice(temperatures[0], temperatures[-1] + dy, dy),
#                 slice(co2ppressures[0], co2ppressures[-1] + dx, dx)]
# z = data_P[:, :-1]
# cmap = plt.get_cmap('rainbow')
# levels = MaxNLocator(nbins=10).tick_values(z.min(), z.max())
# cf = plt.contourf(x[:-1, :-1] + dx / 2.,
#                   y[:-1, :-1] + dy / 2.,
#                   z, levels=levels,
#                   cmap=cmap)
# plt.colorbar(cf)
# plt.title('Amount log10(P) [mol]')
# plt.xlabel('ppCO2 [-]')
# plt.ylabel('T [degC]')
# plt.savefig(results_folder + '/' + 'log10P-vs-T-ppCO2-with-contourf.png', bbox_inches='tight')
# plt.close()
