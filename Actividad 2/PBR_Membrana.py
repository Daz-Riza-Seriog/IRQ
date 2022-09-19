# Code made for Sergio Andrés Díaz Ariza
# 10 September 2022
# License MIT
# IRQ: Python Program-Assigment 1


import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import numpy as np
import seaborn as sns

from Parametros import Parametros
from PBR_DROP import sol_PBR, X_PBR, loc_x_p
import timeit

start = timeit.default_timer()
sns.set()

Rate = Parametros.set_Rate()
Reactor = Parametros.set_Reactor()
ini_state = Parametros.set_initial_state_flux_drop()
Est = Parametros.Estequiometria()


# PBR_calc_f function
def PBR_calc_f(W, x, Rate, Reactor):
    # extract state variables
    F_A = x[0]  # molar flow rate of isopropanol
    F_B = x[1]  # molar flow rate of acetona
    F_C = x[2]  # molar flow rate of hidrogeno
    P = x[3]  # total pressure

    # compute total molar flow rate
    F_tot = F_A + F_B + F_C + Reactor.F_I0

    # compute partial pressures in gas phase
    p_A = (Reactor.P_0 * F_A) / F_tot  # atm
    p_B = (Reactor.P_0 * F_B) / F_tot  # atm
    p_C = (Reactor.P_0 * F_C) / F_tot  # atm

    # compute volumetric flow rate in m^3/s
    R = 8.31446261815324  # Pas*m^3/mol*K
    T = Reactor.T_0  # isothermal operation

    vflow = (F_tot / Reactor.F_tot0) * (P / Reactor.P_0) * (Reactor.T_0 / T) * Reactor.vflow_0  # volumetric flow rate

    # compute surface reaction rate
    r_R = Rate.k_s * (F_A / vflow)

    # Difusión de salida H2
    R_C = Reactor.K * Reactor.a_difu * F_C / vflow

    # Next, evaluate derivative functions.

    # mole balance on A
    dFA_dW = -r_R
    # mole balance on B
    dFB_dW = r_R
    # mole balance on C
    dFC_dW = r_R - R_C

    # ACA DEBEMOS METER LOS DATOS DE PRESION EN TERMINOS DE VOLUMEN #
    # LOS PARAMETROS ESTAN ARRIBA EN LA CLASE PARAMETROS, ES CAMBIARLOS#
    # pressure drop across packed bed
    var_1 = Reactor.beta_0 / Reactor.A_c  # Poner en terminos de volumen
    dP_dV = (-var_1 * (Reactor.P_0 / P) * (F_tot / Reactor.F_tot0))  # Poner en terminos de volumen

    f = [dFA_dW, dFB_dW, dFC_dW, dP_dV]

    return f


t_eval = np.linspace(0, Reactor.V_tot_su, 100000)
sol_mem = solve_ivp(PBR_calc_f, [0, Reactor.V_tot_su], ini_state.x_0, args=(Rate, Reactor), t_eval=t_eval)


def Conv(FA, F_A0):
    X = 1 - (FA / F_A0)
    return X


X_PBR_mem = [Conv(x, Reactor.F_A0) for x in sol_mem.y[0]]


# Find the exact Volume for the Conversion
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


val_x, loc_x = find_nearest(X_PBR_mem, Reactor.X)

# Da = rA0*V/FA0
Da = (Rate.k_s*Reactor.C_A0*sol_mem.t[loc_x])/Reactor.F_A0

print("Volume for PBR-Membrana:", sol_mem.t[loc_x], "[m^3]")
print("Damkohler PBR-Membrana:",Da)

class plot_results:
    fig, axs = plt.subplots(2, 2)

    axs[0, 0].set_title("Conversion vs Volumen")
    axs[0, 0].plot(X_PBR_mem, sol_mem.t, label=r"with Membrane")
    axs[0, 0].plot(X_PBR, sol_PBR.t, label=r"with out Membrane")
    axs[0, 0].set_xlabel('$X$', labelpad=15, fontsize=13)
    axs[0, 0].set_ylabel('Volumen $m^3$', labelpad=8, fontsize=12)
    axs[0, 0].legend()

    axs[0, 1].set_title("Volumen de Rector vs Flujo molar especie $F_{isopropanol}$")
    axs[0, 1].plot(sol_mem.t, sol_mem.y[0], label=r"with Membrane")
    axs[0, 1].plot(sol_PBR.t, sol_PBR.y[0], label=r"with out Membrane")
    axs[0, 1].set_xlabel('Volumen $m^3$', labelpad=15, fontsize=13)
    axs[0, 1].set_ylabel("Flujo molar $[mol/s]$", labelpad=8, fontsize=12)
    axs[0, 1].legend()

    axs[1, 0].set_title("Volumen de Rector vs Flujo molar especie $F_{acetona}$")
    axs[1, 0].plot(sol_mem.t, sol_mem.y[1], label=r"with Membrane")
    axs[1, 0].plot(sol_PBR.t, sol_PBR.y[1], label=r"with out Membrane")
    axs[1, 0].set_xlabel('Volumen $m^3$', labelpad=15, fontsize=13)
    axs[1, 0].set_ylabel("Flujo molar $[mol/s]$", labelpad=8, fontsize=12)
    axs[1, 0].legend()

    axs[1, 1].set_title("Volumen de Rector vs Flujo molar especie $F_{hidrogeno}$")
    axs[1, 1].plot(sol_mem.t, sol_mem.y[2], label=r"with Membrane")
    axs[1, 1].plot(sol_PBR.t, sol_PBR.y[2], label=r"with out Membrane")
    axs[1, 1].set_xlabel('Volumen $m^3$', labelpad=15, fontsize=13)
    axs[1, 1].set_ylabel("Flujo molar $[mol/s]$", labelpad=8, fontsize=12)
    axs[1, 1].legend()
    plt.suptitle("Reactor PBR-Membrane")
    plt.tight_layout()

    plt.figure(2)
    plt.title("PBR-Membrane\nVolumen de Reactor vs Presion ")
    plt.plot(sol_mem.t, sol_mem.y[3], label=r"with Membrane")
    plt.plot(sol_PBR.t, sol_PBR.y[3], label=r"with out Membrane")
    plt.vlines(sol_PBR.t[loc_x_p], sol_mem.y[3][-1], sol_mem.y[3][0], label=r"$X$=0.9-PBR", linestyle="dotted", colors="c")
    plt.vlines(sol_mem.t[loc_x], sol_mem.y[3][-1], sol_mem.y[3][0], label=r"$X$=0.9-PBR-Mem", linestyle="dotted", colors="gold")
    plt.xlabel('Volumen $[m^3]$', labelpad=15, fontsize=13)
    plt.legend()
    plt.ylabel('Presion $[Pas]$', labelpad=8, fontsize=12)
    plt.tight_layout()

    plt.figure(3)
    plt.title("PBR-Membrane\nConversion vs Presion ")
    plt.plot(X_PBR_mem, sol_mem.y[3], label=r"with Membrane")
    plt.plot(X_PBR, sol_PBR.y[3], label=r"with out Membrane")
    plt.xlabel('X$', labelpad=15, fontsize=13)
    plt.legend()
    plt.ylabel('Presion $[Pas]$', labelpad=8, fontsize=12)
    plt.tight_layout()

    plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)
