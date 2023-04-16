from gdtk.gas import GasState, GasModel
"""
gm_1 = GasModel('nitrogen-2sp.lua')
gs_1 = GasState(gm_1)
gs_1.rho = 0.070192050
gs_1.u = 10145002.916967910
gs_1.massf = [0.8692812427881279, 0.13071875721187218]

gs_1.update_thermo_from_rhou()
print("T:", gs_1.T, "p:", gs_1.p)
"""
gm_2 = GasModel('nitrogen-2sp.lua')
gs_2 = GasState(gm_2)
gs_2.T = 6177.351
gs_2.p = 1.455e5
gs_2.massf = [8.692812671004e-01, 1.307187328996e-01]
gs_2.update_thermo_from_pT()
print("rho:", gs_2.rho, "u:", gs_2.u)

gm_3 = GasModel('nitrogen-2sp.lua')
gs_3 = GasState(gm_3)
gs_3.rho = 0.07018354320899349
gs_3.T = 6177.351
gs_3.massf = [8.692812671004e-01, 1.307187328996e-01]
gs_3.update_thermo_from_rhoT()

print("T:", gs_3.T, "p:", gs_3.p, "u:", gs_3.u)
