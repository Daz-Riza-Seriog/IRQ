# Code made for Sergio Andrés Díaz Ariza
# 10 September 2022
# License MIT
# IRQ: Python Program-Assigment 1


import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import numpy as np
import seaborn as sns

from Parametros import Parametros
import timeit

start = timeit.default_timer()
sns.set()

Rate = Parametros.set_Rate()
Reactor = Parametros.set_Reactor()
ini_state = Parametros.set_initial_state_flux_drop()
Est = Parametros.Estequiometria()

# PBR_calc_f function
v_PFR = []


def PBR_calc_f(W, x, Rate, Reactor, v_PFR):
    # extract state variables
    F_A = x[0]  # molar flow rate of isopropanol
    F_B = x[1]  # molar flow rate of acetona
    F_C = x[2]  # molar flow rate of hidrogeno
    # P = x[3]  # total pressure

    # compute total molar flow rate
    F_tot = F_A + F_B + F_C + Reactor.F_I0

    # compute partial pressures in gas phase
    p_A = (Reactor.P_0 * F_A) / F_tot  # atm
    p_B = (Reactor.P_0 * F_B) / F_tot  # atm
    p_C = (Reactor.P_0 * F_C) / F_tot  # atm

    # compute volumetric flow rate in m^3/s
    R = 8.31446261815324  # Pas*m^3/mol*K
    T = Reactor.T_0  # isothermal operation

    vflow = (F_tot / Reactor.F_tot0) * (Reactor.P_0 / Reactor.P_0) * (
            Reactor.T_0 / T) * Reactor.vflow_0  # volumetric flow rate
    v_PFR.append(vflow)

    # compute surface reaction rate
    r_R = Rate.k_s * (F_A / vflow)

    # Next, evaluate derivative functions.

    # mole balance on A
    dFA_dW = -r_R
    # mole balance on B
    dFB_dW = r_R
    # mole balance on C
    dFC_dW = r_R

    # ACA DEBEMOS METER LOS DATOS DE PRESION EN TERMINOS DE VOLUMEN #
    # LOS PARAMETROS ESTAN ARRIBA EN LA CLASE PARAMETROS, ES CAMBIARLOS#
    # pressure drop across packed bed
    var_1 = Reactor.beta_0 / Reactor.A_c  # Poner en terminos de volumen
    dP_dV = (-var_1 * (Reactor.P_0 / Reactor.P_0) * (F_tot / Reactor.F_tot0))  # Poner en terminos de volumen

    f = [dFA_dW, dFB_dW, dFC_dW, dP_dV]

    return f


t_eval = np.linspace(0, Reactor.V_tot_su, 16482)
sol_PFR = solve_ivp(PBR_calc_f, [0, Reactor.V_tot_su], ini_state.x_0, args=(Rate, Reactor, v_PFR), t_eval=t_eval)


def Conv(FA, F_A0):
    X = 1 - (FA / F_A0)
    return X


X_PFR = [Conv(x, Reactor.F_A0) for x in sol_PFR.y[0]]

# Find the exact value for Volume
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


val_x, loc_x = find_nearest(X_PFR, Reactor.X)
Da = (Rate.k_s*Reactor.C_A0*sol_PFR.t[loc_x])/Reactor.F_A0

print("Volume for PBR-Sin Caida de Presion:", sol_PFR.t[loc_x], "[m^3]")
print("Damkohler PBR-Sin caida de Presion:",Da)


class plot_results:
    plt.figure(1)
    plt.title("PBR\nConversion vs Volumen")
    plt.plot(X_PFR, sol_PFR.t)
    plt.xlabel('$X$', labelpad=15, fontsize=13)
    plt.ylabel('Volumen $m^3$', labelpad=8, fontsize=12)
    plt.tight_layout()

    plt.figure(2)
    plt.title("PFR\nVolumen de Reactor vs Flujo molar especie $i$")
    plt.plot(sol_PFR.t, sol_PFR.y[0], label="isopropanol")
    plt.plot(sol_PFR.t, sol_PFR.y[1], label="acetona")
    plt.plot(sol_PFR.t, sol_PFR.y[2], label="hidrogeno")
    plt.xlabel('Volumen $[m^3]$', labelpad=15, fontsize=13)
    plt.legend()
    plt.ylabel('Flujo molar $F_i$ $[mol/s]$', labelpad=8, fontsize=12)
    plt.tight_layout()

    plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)
