# Code made by Sergio Andrés Díaz Ariza
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

class ODE_PBR:

    def __init__(self,x,Rate,Reactor) -> None:
        self.x = x
        self.Rate = Rate
        self.Reactor = Reactor

    # PBR_calc_f function
    def dFdW(self,W,x):
        # extract state variables
        F_A = x[0]  # molar flow rate of CO2
        F_B = x[1]  # molar flow rate of H2
        F_C = x[2]  # molar flow rate of Metanol
        F_D = x[3]  # molar flow rate of Agua
        F_E = x[4]  # molar flow rate of CO
        F_F = x[5]  # molar flow rate of CH4
        r1C = x[6]
        r2A = x[7]
        r3F = x[8]
        P = x[9]  # total pressure
        
        # compute total molar flow rate
        F_tot = F_A + F_B + F_C + F_D + F_E + F_F
        # compute volumetric flow rate in m^3/s
        R = Reactor.R  # Bar*L/mol*K
        T = Reactor.T_0  # isothermal operation
        Z = 1  # Acentric Factor
        
        vflow = (F_tot / Reactor.F_tot0) * (P / Reactor.P_0) * (Reactor.T_0 / T) * Reactor.vflow_0  # volumetric flow rate
        
        # compute partial pressures in gas phase
        P_A = (F_A * Z * R * T) / vflow  # bar
        P_B = (F_B * Z * R * T) / vflow  # bar
        P_C = (F_C * Z * R * T) / vflow  # bar
        P_D = (F_D * Z * R * T) / vflow  # bar
        P_E = (F_E * Z * R * T) / vflow  # bar
        P_F = (F_F * Z * R * T) / vflow  # bar
        
        # pressure drop across packed bed
        var_1 = Reactor.beta_0 / (Reactor.area * (1 - Reactor.f) * Reactor.rc)  # Poner en terminos de volumen
        dP_dW = (var_1 * (Reactor.P_0 / P) * (F_tot / Reactor.F_tot0))  # Poner en terminos de volumen
        
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
        rA = -r1A - r2A - r3A
        rB = -r1B - r2B - r3B
        rC = r1C
        rD = r1D + r2D + r3D
        rE = r2E
        rF = r3F

        # mole balance on CO2
        dFA_dW = rA
        
        # mole balance on H2
        dFB_dW = rB
        # mole balance on C3OH
        dFC_dW = rC
        # mole balance on H20
        dFD_dW = rD
        # mole balance on CO
        dFE_dW = rE
        # mole balance on CH4
        dFF_dW = rF
        
        f = [dFA_dW, dFB_dW, dFC_dW,  dFD_dW,dFE_dW, dFF_dW, r1C, r2A, r3F, dP_dW]

        return f

ODE_1 = ODE_PBR(ini_state.x_0,Rate,Reactor)
#Note we can now simulate two PBR with different condittion without alterate the code and it save the parameters for each one
ODE_2 = ODE_PBR(ini_state.x_1,Rate,Reactor)

t_eval = np.linspace(0, Reactor.V_tot_su, 10000)
sol_PBR = solve_ivp(ODE_2.dFdW,[0, Reactor.V_tot_su], ini_state.x_0, t_eval=t_eval, rtol=1e-12,
                    atol=1e-12)
print(sol_PBR)

# Selectividad r1/r2+r3
S_1 = np.divide(sol_PBR.y[6], np.add(sol_PBR.y[7], sol_PBR.y[8]))
# Selectividad r1/r2
S_2 = np.divide(sol_PBR.y[6], sol_PBR.y[7])
# Selectividad r1/r3
S_3 = np.divide(sol_PBR.y[6], sol_PBR.y[8])

#Rendimiento 1
Y_1_1 = np.divide(sol_PBR.y[2],sol_PBR.y[0][0]-sol_PBR.y[0])

# Find the position of maximum selectivity 1
S_max = np.asarray(S_1[100:]).argmax()
Sel_max = round(S_1[S_max], 4)
Cat_Sel_max = round(sol_PBR.t[S_max], 4)
print("maximo selectividad ´maximo grafica´", Sel_max)

# Find the position of maximum selectivity 2
S_max2 = np.asarray(S_2[100:]).argmax()
Sel_max2 = round(S_2[S_max2], 4)
Cat_Sel_max2 = round(sol_PBR.t[S_max2], 4)
print("maximo selectividad ´maximo grafica´", Sel_max2)

#Find the position of maximum selectivity 1
S_max3 = np.asarray(S_3[100:]).argmax()
Sel_max3=round(S_3[S_max3], 4)
Cat_Sel_max3=round(sol_PBR.t[S_max3], 4)
print("maximo selectividad ´maximo grafica´",Sel_max3)

#Find the position of maximum Global Yield rC3OH_1/rCo2_1
Y_max_1 = np.asarray(Y_1_1[10:]).argmax()
Yil_max_1=round(Y_1_1[Y_max_1], 4)
Time_Yil_max_1=round(sol_PBR.t[Y_max_1], 4)
print("maximo selectividad ´maximo grafica´",Yil_max_1)

# Selectividad from paper equation stop in the value of max selectivity find before
Sel_metanol = sol_PBR.y[2][S_max] / (sol_PBR.y[0][0] - sol_PBR.y[0][S_max])
# Conversion from paper equation  stop in the value of max selectivity find before
Conv_CO2 = (sol_PBR.y[0][0] - sol_PBR.y[0][S_max]) / sol_PBR.y[0][0]
# Production of Metanol that we reach
Prod_Metanol = (sol_PBR.y[2][S_max]) * (32.04 / 1000)

print("Selectividad Metanol ´paper´", Sel_metanol)
print("Conversion CO2", Conv_CO2)
print("Peso de Catalizador usado", sol_PBR.t[S_max], "[Kg]")
print("Produccion Metanol por lote", Prod_Metanol, "Kg/seg")


class plot_results:
    plt.figure(1)
    plt.title("PBR\n$Kg_{cat}$ vs Selectividad")
    plt.plot(sol_PBR.t, S_1, label=r"Selectividad $r_{CH_3OH}/r_{CO}+r_{CH_4}$")
    plt.vlines(sol_PBR.t[S_max], 0, S_1[S_max], label=f"Selectividad Maxima {Sel_max},\nCatalizador {Cat_Sel_max}",
               linestyle="dotted", colors="c")
    plt.xlabel("Catalizador $[Kg]$", labelpad=15, fontsize=13)
    plt.ylabel(r'$r_{selectividad}$', labelpad=8, fontsize=12)
    plt.legend()
    plt.tight_layout()

    plt.figure(2)
    plt.title("PBR\n$Kg_{cat}$ vs Selectividad")
    plt.plot(sol_PBR.t, S_2, label=r"Selectividad $r_{CH_3OH}/r_{CO}$")
    plt.vlines(sol_PBR.t[S_max2], 0, S_1[S_max2], label=f"Selectividad Maxima {Sel_max2},\nCatalizador {Cat_Sel_max2}",
               linestyle="dotted", colors="c")
    plt.xlabel("Catalizador $[Kg]$", labelpad=15, fontsize=13)
    plt.ylabel(r'$r_{selectividad}$', labelpad=8, fontsize=12)
    plt.legend()
    plt.tight_layout()

    plt.figure(3)
    plt.title("PBR\n$Kg_{cat}$ vs Selectividad")
    plt.plot(sol_PBR.t, S_3, label=r"Selectividad $r_{CH_3OH}/r_{CH_4}$")
    plt.vlines(sol_PBR.t[S_max3], 0, S_3[S_max3], label=f"Selectividad Maxima {Sel_max},\nCatalizador {Cat_Sel_max}",
               linestyle="dotted", colors="c")
    plt.xlabel("Catalizador $[Kg]$", labelpad=15, fontsize=13)
    plt.ylabel(r'$r_{selectividad}$', labelpad=8, fontsize=12)
    plt.legend()
    plt.tight_layout()

    plt.figure(4)
    plt.title("PFR\nTiempo vs Rendimiento")
    plt.plot(sol_PBR.t, Y_1_1, label=r"Rendimiento $F_{CH_3OH}/F_{CO_0}-F_{CO}}$")
    plt.vlines(sol_PBR.t[Y_max_1], 0, Y_1_1[Y_max_1], label=f"Rendimiento Maximo {Yil_max_1},\nCatalizador{Time_Yil_max_1}",
               linestyle="dotted", colors="c")
    plt.xlabel("Catalizador $[Kg]$", labelpad=15, fontsize=13)
    plt.ylabel(r'${Rendimiento}$', labelpad=8, fontsize=12)
    plt.legend()
    plt.tight_layout()

    plt.figure(5)
    plt.title("PBR\nCatalizador vs Flujo")
    plt.plot(sol_PBR.t, sol_PBR.y[0], label=r"$F_{CO_2}$")
    plt.plot(sol_PBR.t, sol_PBR.y[1], label=r"$F_{H_2}$")
    plt.plot(sol_PBR.t, sol_PBR.y[2], label=r"$F_{CH_3OH}$")
    plt.plot(sol_PBR.t, sol_PBR.y[3], label=r"$F_{H_2O}$")
    plt.plot(sol_PBR.t, sol_PBR.y[4], label=r"$F_{CO_2}$")
    plt.plot(sol_PBR.t, sol_PBR.y[5], label=r"$F_{CH_4}$")
    plt.xlabel('Catalizador $[Kg]$', labelpad=15, fontsize=13)
    plt.ylabel('Flujo molar $F_i$ $[mol/s]$', labelpad=8, fontsize=12)
    plt.legend()
    plt.tight_layout()

    plt.figure(6)
    plt.plot(sol_PBR.t, sol_PBR.y[6])
    plt.suptitle("PBR\n Catalizador vs $r_{CH_3OH}$")
    plt.ylabel('r $[mol/s*bar^2*Kg_{cat}]$', labelpad=8, fontsize=12)
    plt.xlabel('Catalizador $[Kg]$', labelpad=15, fontsize=13)

    plt.figure(7)
    plt.plot(sol_PBR.t, sol_PBR.y[7])
    plt.suptitle("PBR\n Catalizador vs $r_{CO_2}$")
    plt.ylabel('r $[mol/s*bar^1.5k*Kg_{cat}]$', labelpad=8, fontsize=12)
    plt.xlabel('Catalizador $[Kg]$', labelpad=15, fontsize=13)

    plt.figure(8)
    plt.plot(sol_PBR.t, sol_PBR.y[8])
    plt.suptitle("PBR\n Catalizador vs $r_{CH_4}$")
    plt.xlabel('Catalizador $[Kg]$', labelpad=15, fontsize=13)
    plt.ylabel('r $[mol/s*bar^1*Kg_{cat}]$', labelpad=8, fontsize=12)

    plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)
