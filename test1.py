from gdtk.gas import GasState, GasModel
from copy import deepcopy
import pyeq
import numpy as np
from Algorithms.StoichiometricRatioSolver.molecule import Molecule
from Algorithms.StoichiometricRatioSolver.bulk_reaction import BulkReaction
from Algorithms.DT_1D_V4.models.cell_methods.goal_massf_solver \
    import form_goal_massf_solver_matrices, goal_massf_solver
"""
import plotly.express as px

df = px.data.iris()
fig = px.scatter_3d(df, x='sepal_length', y='sepal_width', z='petal_width',
              color='petal_length', size='petal_length', size_max=18,
              symbol='species', opacity=0.7)

# tight layout
fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
fig.write_html("test.html")
"""
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
"""
"""
gm = GasModel("cea-H2-O2-O-OH-H2O-H-gas-model.lua")
gs = GasState(gm)

gs.p = 33046572
gs.T = 983.15

gs.update_thermo_from_pT()

print(gs.ceaSavedData)
"""
"""
H2 = Molecule(chemical_form = "H2")
O2 = Molecule(chemical_form = "O2")
C = Molecule(chemical_form = "C")
O = Molecule(chemical_form = "O")
CO2 = Molecule(chemical_form = "CO2")
H2O = Molecule(chemical_form = "H2O")

reaction = BulkReaction(species = [H2, O2, C, O, CO2, H2O])

current_massf = [0.0, 0.8, 0.2, 0.0, 0.0, 0.0]
species_names = ["H2", "O2", "C", "O", "CO2", "H2O"]
molar_masses = [2.016, 31.998, 12.011, 15.999, 44.01, 18.01528]
massf_given = {
    "O2"    : 0.0,
    "C"     : 0.0,
    "H2O"    : 0.0
}

LHS_A, RHS_A, source = form_goal_massf_solver_matrices(massf_current = current_massf, \
                    bulk_reaction_parameters = [0.1, massf_given, reaction], \
                    molar_masses = molar_masses, \
                    species_names = species_names)

goal_massf = goal_massf_solver(massf_current = current_massf, \
                    bulk_reaction_parameters = [0.1, massf_given, reaction], \
                    species_names = species_names, \
                    LHS_A = LHS_A, RHS_A = RHS_A, source = source)
print("Carbon test goals massf:", goal_massf)
"""
"""
H2 = Molecule(chemical_form = "H2")
O2 = Molecule(chemical_form = "O2")
OH = Molecule(chemical_form = "OH")
H = Molecule(chemical_form = "H")
O = Molecule(chemical_form = "O")
H2O = Molecule(chemical_form = "H2O")

species_names = ["H2", "O2", "OH", "O", "H", "H2O"]

current_massf = [0.1, 0.9, 0.0, 0.0, 0.0, 0.0]

molar_masses = [2.01588000e-03, 3.19988000e-02, 1.70073400e-02, 1.59994000e-02, 1.00794000e-03, 1.80152800e-02]

massf_given = {
    "O2"     : 0.9,
    "H2"    : 0.1,
    "OH"   : 0.0
}
reaction = BulkReaction(species = [H2, O2, OH, O, H, H2O])

LHS_A, RHS_A, source = form_goal_massf_solver_matrices(massf_current = current_massf, \
                    bulk_reaction_parameters = [0.1, massf_given, reaction], \
                    molar_masses = molar_masses, \
                    species_names = species_names)
print(LHS_A)
print(RHS_A)
print(source)

goal_massf = goal_massf_solver(massf_current = current_massf, \
                    bulk_reaction_parameters = [0.1, massf_given, reaction], \
                    species_names = species_names, \
                    LHS_A = LHS_A, RHS_A = RHS_A, source = source)

print(goal_massf)

"""
###### CEQ testing
"""
mixture = pyeq.EqCalculator(spnames = ["H2", "O2", "OH", "O", "H", "H2O"])
x0 = mixture.YtoX(Y = np.array([0.1, 0.9, 0.0, 0.0, 0.0, 0.0]))
equil_Y = mixture.pt(p = 1e6, T = 3000.0, Xs0 = x0)
print("pyeq equilibrium mass fraction:", mixture.XtoY(X = equil_Y))
mol_masses = mixture.M.tolist()
print("molar masses:", mol_masses)


H2 = Molecule(chemical_form = "H2")
O2 = Molecule(chemical_form = "O2")
OH = Molecule(chemical_form = "OH")
H2O = Molecule(chemical_form = "H2O")
H = Molecule(chemical_form = "H")
O = Molecule(chemical_form = "O")

bulk_reaction = BulkReaction(species = [H2, O2, OH, H2O, H, O])
spcs_names = ["H2", "O2", "OH", "O", "H", "H2O"]
"""
massf_given = {
    "O2"        : 1.13056673e-01,
    "H2"        : 5.08994151e-03,
    "OH"        : 6.45912390e-02,
    "H2O"       : 8.08622451e-01
}
"""
massf_given = {
    "O"        : 0.0,
    "H"        : 0.0,
    "OH"        : 0.0,
    "H2O"       : 1.0
}

current_massf = [0.5, 0.5, 0.0, 0.0, 0.0, 0.0]

LHS_A, RHS_A, source = form_goal_massf_solver_matrices(massf_current = current_massf, \
                    bulk_reaction_parameters = [1.0, massf_given, bulk_reaction], \
                    molar_masses = mol_masses, \
                    species_names = spcs_names)

goal_massf = goal_massf_solver(massf_current = current_massf, \
                    bulk_reaction_parameters = [1.0, massf_given, bulk_reaction], \
                    species_names = spcs_names, \
                    LHS_A = LHS_A, RHS_A = RHS_A, source = source)

print("My solver mass fraciton:", goal_massf)
x = 1.0
for ind, val in enumerate(goal_massf):
    if current_massf[ind] + x * (goal_massf[ind] - current_massf[ind]) < 0.0:
        x = min(x, current_massf[ind] / (current_massf[ind] - goal_massf[ind]))
print("x:", x)
massf_new = (np.array(current_massf) + x * (np.array(goal_massf) - np.array(current_massf))).tolist()
print(massf_new)
"""
