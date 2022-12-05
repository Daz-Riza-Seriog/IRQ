# Code made by Sergio Andrés Díaz Ariza
# 04 November 2022
# License MIT
# IRQ: Python Program-Assigment 5


import matplotlib.pyplot as plt
from Parametros import Parametros
from scipy.integrate import solve_ivp
import numpy as np
import seaborn as sns

import timeit

start = timeit.default_timer()
sns.set()

Reactor = Parametros.set_Reactor()
Rate = Parametros.set_Rate()
ini_state = Parametros.set_initial_state_flux_drop


def PBR_calc_f(W, x, Rate, Reactor):
    # extract state variables
    FA = x[0]  # Flujo molar de CO2
    FB = x[1]  # Flujo molar de H2
    FC = x[2]  # Flujo molar de CH4
    FD = x[3]  # Flujo molar de H2O
    FE = x[4]  # Flujo molar de CO
    P = x[5]  # Presion
    Ta = x[6]
    T = x[7]
    rA1 = x[8]
    rA2 = x[9]
    rA3 = x[10]

    Tref = 623  # Temperatura de Ref for Absorbance [K]

    # compute total molar flow rate
    FT = FA + FB + FC + FD + FE  # mol/s

    # volumetric flow
    vflow = (FT / Reactor.FT0) * (P / Reactor.P0) * (Reactor.T0 / T) * Reactor.vflow_0  # volumetric flow rate

    PCO2 = (FA * Reactor.Rg * T) / vflow  # Presiones parciales CO2
    PH2 = (FB * Reactor.Rg * T) / vflow  # Presiones parciales H2
    PCH4 = (FC * Reactor.Rg * T) / vflow  # Presiones parciales CH4
    PH2O = (FD * Reactor.Rg * T) / vflow  # Presiones parciales H2O
    PCO = (FE * Reactor.Rg * T) / vflow  # Presiones parciales CO

    # Entalpia de Rxn Varia
    H_rxn_1 = Rate.H_rxn_1_est + Rate.D_Cp_1 * (T - Rate.T_est)
    H_rxn_2 = Rate.H_rxn_2_est + Rate.D_Cp_2 * (T - Rate.T_est)
    H_rxn_3 = Rate.H_rxn_3_est + Rate.D_Cp_3 * (T - Rate.T_est)

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

    # pressure drop across packed bed
    flujo_m_tot = (
                          FA * Reactor.PM_A + FB * Reactor.PM_B + FC * Reactor.PM_C + FD * Reactor.PM_D + FE * Reactor.PM_E) / 1000  # [kg/s]
    G = flujo_m_tot / Reactor.Ac  # [kg/m^2 * s]
    ter1 = (150 * (1 - Reactor.f) * Reactor.mu) / Reactor.Dp
    ter2 = (G * (1 - Reactor.f)) / (Reactor.ro * Reactor.Dp * (Reactor.f ** 3))
    beta_0 = ter2 * (ter1 + 1.75 * G)
    var_1 = -beta_0 / (Reactor.Ac * (1 - Reactor.f) * Reactor.rc)

    dFAdW = -rA1 - rA2
    dFBdW = -(4 * rA1) - rA2 - (3 * rA3)
    dFCdW = rA1 + rA3
    dFDdW = (2 * rA1) + rA2 + rA3
    dFEdW = rA2 - rA3
    dPdW = (var_1 * (Reactor.P0 / P) * (FT / Reactor.FT0))  # Poner en terminos de volumen

    # Co-corriente
    Qgen = -rA1 * H_rxn_1 + rA2 * H_rxn_2 - rA3 * H_rxn_3
    Qrem = (Rate.Ua / Reactor.rc) * (T - Ta)  # TODO corregir por que es rho_b
    Cpmix = (FA * Rate.Cp_CO2) + (FB * Rate.Cp_H2) + (FC * Rate.Cp_CH4) + (FD * Rate.Cp_H2O) + (FE * Rate.Cp_CO)
    dTadV = -(Rate.Ua / Reactor.rc) * (T - Ta) / (Reactor.F_cool * Rate.Cp_cool)  # Co-Current Balance
    dTdW = (Qgen - Qrem) / Cpmix

    f = [dFAdW, dFBdW, dFCdW, dFDdW, dFEdW, dPdW, dTadV, dTdW, rA1, rA2,rA3]

    return f


t_eval = np.linspace(0, Reactor.Wf, 1000)
sol_PBR = solve_ivp(PBR_calc_f, [0, Reactor.Wf], ini_state.x_0, args=(Rate, Reactor), t_eval=t_eval, rtol=1e-12,
                    atol=1e-12)


class plot_results:

    fig, axs = plt.subplots(2, 2)

    axs[0, 0].set_title("PBR Counter Current\nCatalizador vs Flujo")
    axs[0, 0].plot(sol_PBR.t, sol_PBR.y[0], label=r"$F_{CO_2}$")
    axs[0, 0].plot(sol_PBR.t, sol_PBR.y[1], label=r"$F_{H_2}$")
    axs[0, 0].plot(sol_PBR.t, sol_PBR.y[2], label=r"$F_{CH_4}$")
    axs[0, 0].plot(sol_PBR.t, sol_PBR.y[3], label=r"$F_{H_2O}$")
    axs[0, 0].plot(sol_PBR.t, sol_PBR.y[4], label=r"$F_{CO}$")
    axs[0, 0].set_xlabel('Catalizador $[Kg]$', labelpad=15, fontsize=13)
    axs[0, 0].set_ylabel('Flujo molar $F_i$ $[mol/s]$', labelpad=8, fontsize=12)

    axs[0, 0].legend()

    axs[0, 1].set_title("PBR Counter Current\n Catalizador vs Caida Presion")
    axs[0, 1].plot(sol_PBR.t, sol_PBR.y[5])
    axs[0, 1].set_xlabel('Catalizador $[Kg]$', labelpad=15, fontsize=13)
    axs[0, 1].set_ylabel('P $[Pa]$', labelpad=8, fontsize=12)

    axs[0, 1].legend()

    axs[1, 0].set_title("PBR Counter Current\n Temperatura Ta vs T")
    axs[1, 0].plot(sol_PBR.t, sol_PBR.y[6], label=r"$T_a$")
    axs[1, 0].plot(sol_PBR.t, sol_PBR.y[7], label=r"$T$")
    axs[1, 0].set_xlabel('Catalizador $[Kg]$', labelpad=15, fontsize=13)
    axs[1, 0].set_ylabel('T $[K]$', labelpad=8, fontsize=12)

    axs[1, 0].legend()

    axs[1, 1].set_title("PBR Counter Current\n Temperatura Ta vs T")
    axs[1, 1].plot(sol_PBR.t, sol_PBR.y[8], label=r"$r_1$")
    axs[1, 1].plot(sol_PBR.t, sol_PBR.y[9], label=r"$r_2$")
    axs[1, 1].plot(sol_PBR.t, sol_PBR.y[10], label=r"$r_3$")
    axs[1, 1].set_xlabel('Catalizador $[Kg]$', labelpad=15, fontsize=13)
    axs[1, 1].set_ylabel('r', labelpad=8, fontsize=12)

    axs[1, 1].legend()
    plt.tight_layout()
plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)
