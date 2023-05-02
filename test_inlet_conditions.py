from gdtk.gas import GasModel, GasState
import pyeq
import numpy as np

psia_to_pa_conversion = 6894.76

lb_to_kg_conversion = 0.453592

M_H2 = 2.016 * 1e-3
M_O2 = 31.998 * 1e-3
M_H = 1.00784 * 1e-3
M_O = 15.999 * 1e-3
M_OH = 17.008 * 1e-3
M_H2O = 18.01528 * 1e-3

def fahrenheit_to_kelvin(F):
    return (F - 32.0) * 5.0 / 9.0 + 273.15

####### oxidiser pre-burner state

m_dot_O2_ob = 24.9 #lb/s
m_dot_H2_ob = 41.1 #lb/s

X_H2_ob = (m_dot_H2_ob / M_H2) / ((m_dot_H2_ob / M_H2) + (m_dot_O2_ob / M_O2)) #= 
X_O2_ob = (m_dot_O2_ob / M_O2) / ((m_dot_H2_ob / M_H2) + (m_dot_O2_ob / M_O2)) #=
print("X_H2_ob:", X_H2_ob)
print("X_O2_ob:", X_O2_ob)


mixture_ob = pyeq.EqCalculator(spnames = ["H2", "O2", "H", "O", "OH", "H2O"])

gm_ob = GasModel("H2-O2-6sp-thermally-perfect-gas-model.lua")
gs_ob = GasState(gm_ob)

molef_init_ob = np.array([X_H2_ob, X_O2_ob, 0.0, 0.0, 0.0, 0.0])
#print(np.shape(molef_init))
massf_init_ob = mixture_ob.XtoY(X = molef_init_ob)
#print(np.shape(massf_init))
gs_ob.massf = massf_init_ob.tolist()
gs_ob.p = 3099.0 * psia_to_pa_conversion
gs_ob.T = fahrenheit_to_kelvin(F = 1087)

gs_ob.update_thermo_from_pT()

cea_gm = GasModel("cea-H2-O2-O-OH-H2O-H-gas-model.lua")
cea_gs = GasState(cea_gm)

#cea_gs.p = 3099.0 * psia_to_pa_conversion
#cea_gs.T = fahrenheit_to_kelvin(F = 728.0)
#cea_gs.update_thermo_from_pT()
#print("cea data:", cea_gs.ceaSavedData)

cea_gs.rho = gs_ob.rho
cea_gs.u = gs_ob.u
cea_gs.update_thermo_from_rhou()
print("cea data:", cea_gs.ceaSavedData)

eq_molef_ob = mixture_ob.pt(p = 3099.0 * psia_to_pa_conversion, T = fahrenheit_to_kelvin(F = 728.0), Xs0 = molef_init_ob)
#print(np.shape(eq_molef))
eq_massf_ob = mixture_ob.XtoY(X = eq_molef_ob)
#print(np.shape(eq_massf))

#eq_molef = mixture.rhou(rho = gs.rho, u = gs.u, Xs0 = molef_init)
#print(np.shape(eq_molef))
#eq_massf = mixture.XtoY(X = eq_molef)

print("eq_massf oxidiser pre-burner:", eq_massf_ob)

##### Fuel pre-burner state

m_dot_O2_fb = 66.9 #lb/s
m_dot_H2_fb = 77.7 # lb/s

X_H2_fb = (m_dot_H2_fb / M_H2) / ((m_dot_H2_fb / M_H2) + (m_dot_O2_fb / M_O2))
X_O2_fb = (m_dot_O2_fb / M_O2) / ((m_dot_H2_fb / M_H2) + (m_dot_O2_fb / M_O2))

mixture_fb = pyeq.EqCalculator(spnames = ["H2", "O2", "H", "O", "OH", "H2O"])

molef_init_fb = np.array([X_H2_fb, X_O2_fb, 0.0, 0.0, 0.0, 0.0])

eq_molef_fb = mixture_fb.pt(p = 3091.0 * psia_to_pa_conversion, T = fahrenheit_to_kelvin(F = 1087.0), Xs0 = molef_init_fb)

eq_massf_fb = mixture_fb.XtoY(X = eq_molef_fb)

print("eq_massf fuel pre-burner:", eq_massf_fb)

#### Initial condition state

m_dot_O2_total = m_dot_O2_fb + 839.0 + m_dot_O2_ob #lb/s
m_dot_H2_total = m_dot_H2_fb + m_dot_H2_ob #lb/s

mixture_init = pyeq.EqCalculator(spnames = ["H2", "O2", "H", "O", "OH", "H2O"])

X_H2_init = (m_dot_H2_total / M_H2) / ((m_dot_H2_total / M_H2) + (m_dot_O2_total / M_O2))
X_O2_init = (m_dot_O2_total / M_O2) / ((m_dot_H2_total / M_H2) + (m_dot_O2_total / M_O2))

P_c = 2871.0 #psia
T_c = 6000.0 #deg F

molef_init_cc = np.array([X_H2_init, X_O2_init, 0.0, 0.0, 0.0, 0.0])

eq_molef_init = mixture_fb.pt(p = P_c * psia_to_pa_conversion, T = fahrenheit_to_kelvin(F = T_c), Xs0 = molef_init_cc)

eq_massf_init = mixture_fb.XtoY(X = eq_molef_init)

print("eq_massf combustion chamber:", eq_massf_init)

### mixture inlet condition

