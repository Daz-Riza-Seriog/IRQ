# Code made by Sergio Andrés Díaz Ariza
# 10 September 2022
# License MIT
# IRQ: Python Program-Assigment 1


import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import scipy.integrate as integrate  # Usamos para integrar
import numpy as np
import seaborn as sns

from Parametros import Parametros
import timeit

start = timeit.default_timer()
sns.set()

Rate = Parametros.set_Rate()  # llamamos las K'as de las Reaciones
Reactor = Parametros.set_Reactor()  # llamamos los Parametros del reactor
ini_state = Parametros.set_initial_state_batch()  # llamamos los Estados iniciales del reactor
Est = Parametros.Estequiometria()  # llamamos Realaciones de Estequiometria


# PBR_calc_f function

def BATCH_calc(t, x, Rate, Reactor, Est):
    # extract state variables
    C_A = x[0]  # molar flow rate of CO2
    C_B = x[1]  # molar flow rate of H2
    C_C = x[2]  # molar flow rate of Metanol
    C_D = x[3]  # molar flow rate of Agua
    C_E = x[4]  # molar flow rate of CO
    C_F = x[5]  # molar flow rate of CH4
    r1C = x[6]
    r2A = x[7]
    r3F = x[8]

    # In partial Presures with isotermic Reactor
    P_A = C_A * (Reactor.R * Reactor.T_0)
    P_B = C_B * (Reactor.R * Reactor.T_0)
    P_C = C_C * (Reactor.R * Reactor.T_0)
    P_D = C_D * (Reactor.R * Reactor.T_0)
    P_E = C_E * (Reactor.R * Reactor.T_0)
    P_F = C_F * (Reactor.R * Reactor.T_0)

    # print(P_A)

    # Compute reaction rate
    # Inhibition Parameter
    inh = (1 + (Rate.K_ads_co2 * P_A) + np.sqrt(Rate.K_ads_h2 * P_B)) ** 2

    # Rate law 1
    num_r1 = (P_A * (P_B ** 3) - ((P_C * P_D) / Rate.k_eq_1)) / (P_B ** 2)  # Numerador velocidad de reacción 1
    r1C = (Rate.k_vel_1 * num_r1) / inh  # mol/kg cat*s
    r1A = Est.a_c * r1C  # mol/kg cat*s
    r1B = Est.b_c * r1C  # mol/kg cat*s
    r1D = Est.d_c * r1C  # mol/kg cat*s

    # Rate law 2
    num_r2 = ((P_A * P_B) - ((P_E * P_D) / Rate.k_eq_2)) / (np.sqrt(P_B))  # Numerador velocidad de reacción 2
    r2A = (Rate.k_vel_2 * num_r2) / inh  # mol/kg cat*s
    r2B = (1 / 1) * r2A  # mol/kg cat*s
    r2E = (1 / 1) * r2A  # mol/kg cat*s
    r2D = (1 / 1) * r2A  # mol/kg cat*s

    # Rate law 3
    num_r3 = (np.sqrt(P_A) * np.sqrt(P_B) * (
            1 - ((P_F * (P_D ** 2)) / (P_A * (P_B ** 4) * Rate.k_eq_3))))  # Numerador velocidad de reacción 3

    r3F = (Rate.k_vel_3 * num_r3) / inh  # mol/kg cat*s
    r3A = (1 / 1) * r3F  # mol/kg cat*s
    r3B = (4 / 1) * r3F  # mol/kg cat*s
    r3D = (2 / 1) * r3F  # mol/kg cat*s

    # Net law rate for each component
    rA = -r1A - r3A - r2A
    rB = -r1B - r3B - r2B
    rC = r1C
    rD = r1D + r3D + r2D
    rE = r2E
    rF = r3F

    # mole balance on CO2
    dCA_dt = rA
    # mole balance on H2
    dCB_dt = rB
    # mole balance on C3OH
    dCC_dt = rC
    # mole balance on H20
    dCD_dt = rD
    # mole balance on CO
    dCE_dt = rE
    # mole balance on CH4
    dCF_dt = rF

    f = [dCA_dt, dCB_dt, dCC_dt, dCD_dt, dCE_dt, dCF_dt, r1C, r2A, r3F]

    return f


# Solucionamos las EDOs
t_eval = np.linspace(0, Reactor.t_batch, 1000)
sol = solve_ivp(BATCH_calc, [0, Reactor.t_batch], ini_state.x_0, args=(Rate, Reactor, Est), t_eval=t_eval,rtol=1e-12,atol=1e-12)

# Selectividad r1/r2+r3
S_1 = np.divide(sol.y[6], np.add(sol.y[7],sol.y[8]))
#Rendimiento 1
Y_1_1 = np.divide(sol.y[2],sol.y[0][0]-sol.y[0])

#Find the position of maximum selectivity
S_max = np.asarray(S_1[10:]).argmax()
Sel_max=round(S_1[S_max], 4)
Time_Sel_max=round(sol.t[S_max], 4)
print("maximo selectividad ´maximo grafica´",Sel_max)

#Find the position of maximum Global Yield rC3OH_1/rCo2_1
Y_max_1 = np.asarray(Y_1_1[10:]).argmax()
Yil_max_1=round(Y_1_1[Y_max_1], 4)
Time_Yil_max_1=round(sol.t[Y_max_1], 4)
print("maximo selectividad ´maximo grafica´",Yil_max_1)


# Selectividad from paper equation stop in the value of max selectivity find before
Sel_metanol = sol.y[2][S_max] / (sol.y[0][0] - sol.y[0][S_max])
# Conversion from paper equation  stop in the value of max selectivity find before
Conv_CO2 = (sol.y[0][0] - sol.y[0][S_max]) / sol.y[0][0]
# Volume of Batach at maximum selectivity
Vol_batch = Reactor.F_A0 * sol.t[S_max] / sol.y[0][0]
# Production of Metanol that we reach
Prod_Metanol = (sol.y[2][S_max]) * (32.04 / 1000)

print("Selectividad Metanol ´paper´", Sel_metanol)
print("Conversion CO2", Conv_CO2)
print("Volumen de Batch", Vol_batch, "[m^3]")
print("Produccion Metanol por lote", Prod_Metanol, "Kg/lote {:.4f} seg".format(sol.t[S_max]))



class plot_results:
    plt.figure(1)
    plt.title("BATCH\nTiempo vs Selectividad")
    plt.plot(sol.t, S_1, label=r"Selectividad $r_{CH_3OH}/r_{CO}+r_{CH_4}$")
    plt.vlines(sol.t[S_max], 0, S_1[S_max], label=f"Selectividad Maxima {Sel_max},\nTiempo{Time_Sel_max}", linestyle="dotted", colors="c")
    plt.xlabel("time $[s]$", labelpad=15, fontsize=13)
    plt.ylabel(r'$r_{selectividad}$', labelpad=8, fontsize=12)
    plt.legend()
    plt.tight_layout()

    plt.figure(2)
    plt.title("BATCH\nTiempo vs Rendimiento")
    plt.plot(sol.t, Y_1_1, label=r"Rendimiento $F_{CH_3OH}/F_{CO_0}-F_{CO}}$")
    plt.vlines(sol.t[Y_max_1], 0, Y_1_1[Y_max_1], label=f"Rendimiento Maximo {Yil_max_1},\nTiempo{Time_Yil_max_1}",
               linestyle="dotted", colors="c")
    plt.xlabel("time $[s]$", labelpad=15, fontsize=13)
    plt.ylabel(r'${Rendimiento}$', labelpad=8, fontsize=12)
    plt.legend()
    plt.tight_layout()

    plt.figure(3)
    plt.title("BATCH\nTime vs Concentracion")
    plt.plot(sol.t, sol.y[0], label=r"$CO_2$")
    plt.plot(sol.t, sol.y[1], label=r"$H_2$")
    plt.plot(sol.t, sol.y[2], label=r"$CH_3OH$")
    plt.plot(sol.t, sol.y[3], label=r"$H_2O$")
    plt.plot(sol.t, sol.y[4], label=r"$CO_2$")
    plt.plot(sol.t, sol.y[5], label=r"$CH_4$")
    plt.xlabel('time $[s]$', labelpad=15, fontsize=13)
    plt.ylabel('Concentracion $C_i$ $[mol/m^3]$', labelpad=8, fontsize=12)
    plt.legend()
    plt.tight_layout()

    plt.figure(4)
    plt.plot(sol.t, sol.y[6])
    plt.suptitle("BATCH\n Time vs $r_{CH_3OH}$")
    plt.ylabel('r $[mol/s*bar^2*Kg_{cat}]$', labelpad=8, fontsize=12)
    plt.xlabel('time $[s]$', labelpad=15, fontsize=13)

    plt.figure(5)
    plt.plot(sol.t, sol.y[7])
    plt.suptitle("BATCH\n Time vs $r_{CO_2}$")
    plt.ylabel('r $[mol/s*bar^1.5k*Kg_{cat}]$', labelpad=8, fontsize=12)
    plt.xlabel('time $[s]$', labelpad=15, fontsize=13)

    plt.figure(6)
    plt.plot(sol.t, sol.y[8])
    plt.suptitle("BATCH\n Time vs $r_{CH_4}$")
    plt.xlabel('time $[s]$', labelpad=15, fontsize=13)
    plt.ylabel('r $[mol/s*bar^1*Kg_{cat}]$', labelpad=8, fontsize=12)

    plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)
