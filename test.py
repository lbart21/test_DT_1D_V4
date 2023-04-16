from gdtk.gas import GasModel, GasState
import numpy as np
import matplotlib.pyplot as plt

"""
comb_chamber_init_gm = GasModel("H2-O2-H2O-thermally-perfect-gas-model.lua")
comb_chamber_init_gs = GasState(comb_chamber_init_gm)
comb_chamber_init_molef = [2/3, 1/3, 0.0]
comb_chamber_init_gs.p = 1.9795e7
#comb_chamber_init_gs.T = 3588.706
comb_chamber_init_gs.T = 300.0
print()
"""
n_points = 100
T_list = np.linspace(500, 8000, n_points)
massf_list = []
p = 22587000.0
fig = plt.figure(figsize=(15, 5))
for i in range(n_points):
    gm = GasModel("cea-H2-O2-O-OH-H2O-H-gas-model.lua")
    gs = GasState(gm)
    gs.p = p
    gs.T = T_list[i]
    gs.update_thermo_from_pT()
    massf_list.append(list(gs.ceaSavedData["massf"].values()))
    
transposed_massf = list(map(list, zip(*massf_list)))
#print(transposed_massf)
for j in range(len(transposed_massf)):
    plt.scatter(T_list, transposed_massf[j], label = list(gs.ceaSavedData["massf"].keys())[j])

plt.grid()
plt.legend()
plt.show()


gm_test = GasModel("cea-H2-O2-O-OH-H2O-H-gas-model.lua")
gs_test = GasState(gm_test)
gs_test.p = p
gs_test.T = 300.0
gs_test.update_thermo_from_pT()
print(gs_test.ceaSavedData["massf"])

