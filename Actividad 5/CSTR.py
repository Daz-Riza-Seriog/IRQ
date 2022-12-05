# Code made by Sergio Andrés Díaz Ariza
# 13 November 2022
# License MIT
# IRQ: Python Program-Assigment 5


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from Parametros import Parametros
import warnings

warnings.filterwarnings('ignore', 'The iteration is not making good progress')
from scipy.optimize import fsolve
from scipy.optimize import least_squares

import timeit

start = timeit.default_timer()
sns.set()

Reactor = Parametros.set_Reactor()
Rate = Parametros.set_Rate()
ini_state = Parametros.set_initial_state_flux_cstr


def fun_CSTRT_T(T, F_guess, Rate, Reactor):
    # Function that solves the mole balance equations for a given T and computes the R(T) and G(T)

    Tref = 623  # Temperatura de Ref for Absorbance [K]

    # Constantes Equilibrio
    k_eq_1 = Rate.k_eq_est1 * np.exp(
        (Rate.H_rxn_1_est / Rate.R) * ((1 / Rate.T_est) - (1 / T)))  # Constante de equilibrio 1
    k_eq_2 = Rate.k_eq_est2 * np.exp(
        (Rate.H_rxn_2_est / Rate.R) * ((1 / Rate.T_est) - (1 / T)))  # Constante de equilibrio 2
    k_eq_3 = Rate.k_eq_est3 * np.exp(
        (Rate.H_rxn_3_est / Rate.R) * ((1 / Rate.T_est) - (1 / T)))  # Constante de equilibrio 3

    # # CONSTANTES DE VELOCIDAD PARA CADA REACCION REFERENCIA A T OPERACION
    k_CO2Meth = Rate.kf_1 * np.exp(-Rate.Ea_1 / (Rate.R * T))  # Constante cinética [mol min-1 g-1]
    k_RWGS = Rate.kf_2 * np.exp(-Rate.Ea_2 / (Rate.R * T))  # Constante cinética [mol min-1 g-1]
    k_COMeth = Rate.kf_3 * np.exp(-Rate.Ea_3 / (Rate.R * T))  # Constante cinética [mol min-1 g-1]

    # Constantes de Absorcion
    Kabs_H2 = Rate.KAb_H2 * np.exp(
        (Rate.QAb_H2 / Rate.R) * ((1 / Tref) - (1 / T)))  # bar-1  Constante de absorción obtenida del artículo
    Kabs_CO2 = Rate.KAb_CO2 * np.exp(
        (Rate.QAb_CO2 / Rate.R) * ((1 / Tref) - (1 / T)))  # bar-1  Constante de absorción obtenida del artículo
    Kabs_H2O = Rate.KAb_H2O * np.exp(
        (Rate.QAb_H2O / Rate.R) * ((1 / Tref) - (1 / T)))  # bar-1  Constante de absorción obtenida del artículo
    Kabs_CO = Rate.KAb_CO * np.exp(
        (Rate.QAb_CO / Rate.R) * ((1 / Tref) - (1 / T)))  # bar-1  Constante de absorción obtenida del artículo

    """ Kinetic parameters
    """
    # Entalpia de Rxn Varia
    H_rxn_1 = Rate.H_rxn_1_est + Rate.D_Cp_1 * (T - Rate.T_est)
    H_rxn_2 = Rate.H_rxn_2_est + Rate.D_Cp_2 * (T - Rate.T_est)
    H_rxn_3 = Rate.H_rxn_3_est + Rate.D_Cp_3 * (T - Rate.T_est)

    """ Solution for CA concentrations 
    fsolve(func, x0[, args, fprime, ...])
    Find the roots of a function.
    """
    args = (
        T, Rate, Reactor, k_eq_1, k_eq_2, k_eq_3, H_rxn_1, H_rxn_2, H_rxn_3, k_CO2Meth, k_RWGS, k_COMeth, Kabs_CO2,
        Kabs_H2,
        Kabs_H2O, Kabs_CO)
    #bounds = ((0, 0), (10, 10))
    solution = fsolve(CSTR_main, F_guess, args=args, maxfev=10000, xtol=1e-13)
    #solution = least_squares(CSTR_main, F_guess, bounds=bounds, args=args, xtol=1e-13, gtol=1e-13, ftol=1e-13)

    FA = solution[0]  # Flujo molar de CO2
    FB = solution[1]  # Flujo molar de H2
    FC = solution[2]  # Flujo molar de CH4
    FD = solution[3]  # Flujo molar de H2O
    FE = solution[4]  # Flujo molar de CO

    # compute total molar flow rate
    FT = FA + FB + FC + FD + FE  # mol/s

    # volumetric flow
    vflow = (FT / Reactor.FT0) * (Reactor.P0 / Reactor.P0) * (Reactor.T0 / T) * Reactor.vflow_0  # volumetric flow rate

    PCO2 = (FA * Reactor.Rg * T) / vflow  # Presiones parciales CO2
    PH2 = (FB * Reactor.Rg * T) / vflow  # Presiones parciales H2
    PCH4 = (FC * Reactor.Rg * T) / vflow  # Presiones parciales CH4
    PH2O = (FD * Reactor.Rg * T) / vflow  # Presiones parciales H2O
    PCO = (FE * Reactor.Rg * T) / vflow  # Presiones parciales CO

    # VELOCITY LAWS
    # Velocidad 1
    num_rprima1 = (k_CO2Meth * Kabs_H2 * Kabs_CO2 * PH2 * PCO2) * (
            1 - (PCH4 * (PH2O ** 2) / ((PH2 ** 4) * PCO2 * k_eq_1)))
    den_rprima1 = ((1 + (Kabs_CO2 * PCO2) + (Kabs_H2 * PH2) + (Kabs_H2O * PH2O) + (
            Kabs_CO * PCO)) ** 2)  # velocidad de reacción [mol/g min]
    rprima1 = (num_rprima1 / den_rprima1) * (1000 / 60)  # mol/kg s
    rA1 = rprima1  # mol/kg s

    # Velocidad 2
    num_rprima2 = (k_RWGS * Kabs_CO2 * PCO2) * (1 - (PCO * PH2O / (PH2 * PCO2 * k_eq_2)))
    den_rprima2 = (1 + (Kabs_CO2 * PCO2) + (Kabs_H2 * PH2) + (Kabs_H2O * PH2O) + (
            Kabs_CO * PCO))  # velocidad de reacción [mol/g min]
    rprima2 = (num_rprima2 / den_rprima2) * (1000 / 60)  # mol/kg s
    rA2 = rprima2  # mol/kg s

    # Velocidad 3
    num_rprima3 = (k_COMeth * Kabs_H2 * Kabs_CO * PH2 * PCO) * (1 - (PCH4 * PH2O / ((PH2 ** 3) * PCO * k_eq_3)))
    den_rprima3 = ((1 + (Kabs_CO2 * PCO2) + (Kabs_H2 * PH2) + (Kabs_H2O * PH2O) + (
            Kabs_CO * PCO)) ** 2)  # velocidad de reacción [mol/g min]
    rprima3 = (num_rprima3 / den_rprima3) * (1000 / 60)  # mol/kg s
    rA3 = rprima3  # mol/kg s

    """ Define R(T) term from energy balance
    """

    RT = Rate.Ua * (Reactor.Ta0 - T) - Reactor.FA0 * Rate.Cp_CO2 * (T - Reactor.T0) - Reactor.FB0 * Rate.Cp_H2 * (
            T - Reactor.T0) - Reactor.FE0 * Rate.Cp_CO * (T - Reactor.T0)
    RT = - RT

    """ Define G(T) term from energy balance
    """
    print(FE)

    GT = Reactor.V_cstr * (rA1 * H_rxn_1 + rA2 * H_rxn_2 + rA3 * H_rxn_3)

    # print ("T =", round(T,5),"K""\n"
    #        "CA =", round(CA,5),"\n"
    #        "CB =", round(CB,5),"\n"
    #        "CC =", round(CC,5),"\n"
    #        "RT =", round(RT,5),"\n"
    #        "GT =", round(GT,5),"\n")

    return RT, GT, FA, FB, FC, FD, FE


def CSTR_main(F, *args):
    (
    T, Rate, Reactor, k_eq_1, k_eq_2, k_eq_3, H_rxn_1, H_rxn_2, H_rxn_3, k_CO2Meth, k_RWGS, k_COMeth, Kabs_CO2, Kabs_H2,
    Kabs_H2O, Kabs_CO) = args
    # extract state variables
    FA = F[0]  # Flujo molar de CO2
    FB = F[1]  # Flujo molar de H2
    FC = F[2]  # Flujo molar de CH4
    FD = F[3]  # Flujo molar de H2O
    FE = F[4]  # Flujo molar de CO

    # compute total molar flow rate
    FT = FA + FB + FC + FD + FE  # mol/s

    # volumetric flow
    vflow = (FT / Reactor.FT0) * (Reactor.P0 / Reactor.P0) * (Reactor.T0 / T) * Reactor.vflow_0  # volumetric flow rate

    PCO2 = (FA * Reactor.Rg * T) / vflow  # Presiones parciales CO2
    PH2 = (FB * Reactor.Rg * T) / vflow  # Presiones parciales H2
    PCH4 = (FC * Reactor.Rg * T) / vflow  # Presiones parciales CH4
    PH2O = (FD * Reactor.Rg * T) / vflow  # Presiones parciales H2O
    PCO = (FE * Reactor.Rg * T) / vflow  # Presiones parciales CO

    # VELOCITY LAWS
    # Velocidad 1
    num_rprima1 = (k_CO2Meth * Kabs_H2 * Kabs_CO2 * PH2 * PCO2) * (
            1 - (PCH4 * (PH2O ** 2) / ((PH2 ** 4) * PCO2 * k_eq_1)))
    den_rprima1 = ((1 + (Kabs_CO2 * PCO2) + (Kabs_H2 * PH2) + (Kabs_H2O * PH2O) + (
            Kabs_CO * PCO)) ** 2)  # velocidad de reacción [mol/g min]
    rprima1 = (num_rprima1 / den_rprima1) * (1000 / 60)  # mol/kg s
    rA1 = rprima1  # mol/kg s

    # Velocidad 2
    num_rprima2 = (k_RWGS * Kabs_CO2 * PCO2) * (1 - (PCO * PH2O / (PH2 * PCO2 * k_eq_2)))
    den_rprima2 = (1 + (Kabs_CO2 * PCO2) + (Kabs_H2 * PH2) + (Kabs_H2O * PH2O) + (
            Kabs_CO * PCO))  # velocidad de reacción [mol/g min]
    rprima2 = (num_rprima2 / den_rprima2) * (1000 / 60)  # mol/kg s
    rA2 = rprima2  # mol/kg s

    # Velocidad 3
    num_rprima3 = (k_COMeth * Kabs_H2 * Kabs_CO * PH2 * PCO) * (1 - (PCH4 * PH2O / ((PH2 ** 3) * PCO * k_eq_3)))
    den_rprima3 = ((1 + (Kabs_CO2 * PCO2) + (Kabs_H2 * PH2) + (Kabs_H2O * PH2O) + (
            Kabs_CO * PCO)) ** 2)  # velocidad de reacción [mol/g min]
    rprima3 = (num_rprima3 / den_rprima3) * (1000 / 60)  # mol/kg s
    rA3 = rprima3  # mol/kg s

    dFAdW = -rA1 - rA2
    dFBdW = -(4 * rA1) - rA2 - (3 * rA3)
    dFCdW = rA1 + rA3
    dFDdW = (2 * rA1) + rA2 + rA3
    dFEdW = rA2 - rA3

    f1 = Reactor.FA0 - FA + dFAdW * Reactor.V_cstr
    f2 = Reactor.FB0 - FB + dFBdW * Reactor.V_cstr
    f3 = Reactor.FC0 - FC + dFCdW * Reactor.V_cstr
    f4 = Reactor.FD0 - FD + dFDdW * Reactor.V_cstr
    f5 = Reactor.FE0 - FE + dFEdW * Reactor.V_cstr

    f = [f1,f2,f3,f4,f5]

    return f


# Input Data
# ========================================================================================================================#

""" Operation Conditions
    T  - Temperature (K)
    nT - Number of points for temperature iteration
    v  - Volumetric flow (dm3/min)
"""

nT = 100
T = np.linspace(Reactor.T0, 723, nT)

# param.To = 283
# param.Ta = 57 + 273

F_guess = np.array(ini_state.x_0)

# Solution
# ========================================================================================================================

count = 0
FA = np.zeros(len(T))
FB = np.zeros(len(T))
FC = np.zeros(len(T))
FD = np.zeros(len(T))
FE = np.zeros(len(T))
GT = np.zeros(len(T))
RT = np.zeros(len(T))
for i in T:
    sol = fun_CSTRT_T(T[count], F_guess, Rate, Reactor)

    F_guess = [sol[2], sol[3], sol[4], sol[5], sol[6]]

    RT[count] = sol[0]
    GT[count] = sol[1]
    FA[count] = sol[2]
    FB[count] = sol[3]
    FC[count] = sol[4]
    FD[count] = sol[5]
    FE[count] = sol[6]

    count += 1

# Output Data
# ========================================================================================================================#

""" Fitting Curve """

plt.figure(1)
plt.plot(T, RT, label='RT')
plt.plot(T, GT, label='GT')
plt.xlabel('T (K)', labelpad=8, fontsize=12)
plt.ylabel('RT and GT', labelpad=8, fontsize=12)
plt.grid(alpha=.4, linestyle='--')
plt.legend()

plt.show()
