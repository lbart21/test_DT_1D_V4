from gdtk.gas import GasState, GasModel
from copy import deepcopy

"""
gm = GasModel("cea-nitrogen-2sp-gas-model.lua")
gs = GasState(gm)
gs.p = 1e5
gs.T = 4000.0

gs.update_thermo_from_pT()
"""
"""
rho = 0.084156
u = 3483870.0
print("rho:", rho)
print("u:", u)

gm_actual = GasModel("cea-nitrogen-2sp-gas-model.lua")
gs_actual = GasState(gm_actual)
gs_actual.rho = deepcopy(rho)
gs_actual.u = deepcopy(u)
gs_actual.update_thermo_from_rhou()

#print("gs data:", gs.ceaSavedData)
print("gs_actual data:", gs_actual.ceaSavedData)
"""

gm_dummy = GasModel("H2-O2-N2-9sp-thermally-perfect-gas-model.lua")
print(gm_dummy.species_names)
gs_dummy = GasState(gm_dummy)
gs_dummy.molef = {'H2':0.5, 'O2':0.5, 'N2':0.0, 'H':0.0, 'O':0.0, 'OH':0.0, 'H2O':0.0, 'HO2':0.0, 'H2O2':0.0}
gs_dummy.p = 325.0e3
gs_dummy.T = 4320.0
gs_dummy.update_thermo_from_pT()
print("rho conserved:", gs_dummy.rho, "u conserved:", gs_dummy.u)

gm_goal = GasModel("cea-H2-O2-N2-9sp-gas-model.lua")
gs_goal = GasState(gm_goal)
gs_goal.rho = gs_dummy.rho
gs_goal.u = gs_dummy.u
gs_goal.update_thermo_from_rhou()

massf_equil = gs_goal.ceaSavedData["massf"]

sum = 0.0
for key in massf_equil.keys():
    sum += massf_equil[key]
for key in massf_equil.keys():
    massf_equil[key] /= sum

print("normalised massf_equil", massf_equil)

print(gs_goal.ceaSavedData)