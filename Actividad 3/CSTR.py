# Code made by Sergio Andrés Díaz Ariza
# 10 September 2022
# License MIT
# IRQ: Python Program-Assigment 1


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.optimize import minimize

from Parametros import Parametros
import timeit

start = timeit.default_timer()
sns.set()

Rate = Parametros.set_Rate()  # llamamos las K'as de las Reaciones
Reactor = Parametros.set_Reactor()  # llamamos los Parametros del reactor
ini_state = Parametros.set_initial_state_flux()  # llamamos los Estados iniciales del reactor
Est = Parametros.Estequiometria()  # llamamos Realaciones de Estequiometria


def min_error_Balance(x, *args):
    V = Reactor.V_cstr
    F_A = x[0]  # molar flow rate of CO2
    F_B = x[1]  # molar flow rate of H2
    F_C = x[2]  # molar flow rate of Metanol
    F_D = x[3]  # molar flow rate of Agua
    F_E = x[4]  # molar flow rate of CO
    F_F = x[5]  # molar flow rate of CH4

    # compute total molar flow rate
    F_tot = F_A + F_B + F_C + F_D + F_E + F_F

    # compute volumetric flow rate in m^3/s
    R = Reactor.R  # Bar*L/mol*K
    T = Reactor.T_0  # isothermal operation
    Z = 1  # Acentric Factor

    vflow = (F_tot / Reactor.F_tot0) * (Reactor.P_0 / Reactor.P_0) * (
            Reactor.T_0 / T) * Reactor.vflow_0  # volumetric flow rate

    # compute partial pressures in gas phase
    P_A = (F_A * Z * R * T) / vflow  # bar
    P_B = (F_B * Z * R * T) / vflow  # bar
    P_C = (F_C * Z * R * T) / vflow  # bar
    P_D = (F_D * Z * R * T) / vflow  # bar
    P_E = (F_E * Z * R * T) / vflow  # bar
    P_F = (F_F * Z * R * T) / vflow  # bar

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
    num_r3 = (np.sqrt(P_A * P_B) * (
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

    return ((Reactor.F_A0 - F_A + rA * V) ** 2) + ((Reactor.F_B0 - F_B + rB * V) ** 2) + (
            (Reactor.F_C0 - F_C + rC * V) ** 2) + ((Reactor.F_D0 - F_D + rD * V) ** 2) + (
                   (Reactor.F_E0 - F_E + rE * V) ** 2) + ((Reactor.F_F0 - F_F + rF * V) ** 2)


def constraint_rC(x, *args):
    F_A = x[0]  # molar flow rate of CO2
    F_B = x[1]  # molar flow rate of H2
    F_C = x[2]  # molar flow rate of Metanol
    F_D = x[3]  # molar flow rate of Agua
    F_E = x[4]  # molar flow rate of CO
    F_F = x[5]  # molar flow rate of CH4
    # compute total molar flow rate
    F_tot = F_A + F_B + F_C + F_D + F_E + F_F
    # compute volumetric flow rate in m^3/s
    R = Reactor.R  # Bar*L/mol*K
    T = Reactor.T_0  # isothermal operation
    Z = 1  # Acentric Factor
    vflow = (F_tot / Reactor.F_tot0) * (Reactor.P_0 / Reactor.P_0) * (
            Reactor.T_0 / T) * Reactor.vflow_0  # volumetric flow rate
    # compute partial pressures in gas phase
    P_A = (F_A * Z * R * T) / vflow  # bar
    P_B = (F_B * Z * R * T) / vflow  # bar
    P_C = (F_C * Z * R * T) / vflow  # bar
    P_D = (F_D * Z * R * T) / vflow  # bar
    P_E = (F_E * Z * R * T) / vflow  # bar
    P_F = (F_F * Z * R * T) / vflow  # bar
    # Compute reaction rate
    # Inhibition Parameter
    inh = (1 + (Rate.K_ads_co2 * P_A) + np.sqrt(Rate.K_ads_h2 * P_B)) ** 2
    # Rate law 1
    num_r1 = (P_A * (P_B ** 3) - ((P_C * P_D) / Rate.k_eq_1)) / (P_B ** 2)  # Numerador velocidad de reacción 1
    r1C = (Rate.k_vel_1 * num_r1) / inh  # mol/kg cat*s
    return r1C


def constraint_rCO2(x, *args):
    F_A = x[0]  # molar flow rate of CO2
    F_B = x[1]  # molar flow rate of H2
    F_C = x[2]  # molar flow rate of Metanol
    F_D = x[3]  # molar flow rate of Agua
    F_E = x[4]  # molar flow rate of CO
    F_F = x[5]  # molar flow rate of CH4
    # compute total molar flow rate
    F_tot = F_A + F_B + F_C + F_D + F_E + F_F
    # compute volumetric flow rate in m^3/s
    R = Reactor.R  # Bar*L/mol*K
    T = Reactor.T_0  # isothermal operation
    Z = 1  # Acentric Factor
    vflow = (F_tot / Reactor.F_tot0) * (Reactor.P_0 / Reactor.P_0) * (
            Reactor.T_0 / T) * Reactor.vflow_0  # volumetric flow rate
    # compute partial pressures in gas phase
    P_A = (F_A * Z * R * T) / vflow  # bar
    P_B = (F_B * Z * R * T) / vflow  # bar
    P_C = (F_C * Z * R * T) / vflow  # bar
    P_D = (F_D * Z * R * T) / vflow  # bar
    P_E = (F_E * Z * R * T) / vflow  # bar
    P_F = (F_F * Z * R * T) / vflow  # bar
    # Compute reaction rate
    # Inhibition Parameter
    inh = (1 + (Rate.K_ads_co2 * P_A) + np.sqrt(Rate.K_ads_h2 * P_B)) ** 2
    # Rate law 2
    num_r2 = ((P_A * P_B) - ((P_E * P_D) / Rate.k_eq_2)) / (np.sqrt(P_B))  # Numerador velocidad de reacción 2
    r2A = (Rate.k_vel_2 * num_r2) / inh  # mol/kg cat*s
    return r2A


def constraint_rCH4(x, *args):
    F_A = x[0]  # molar flow rate of CO2
    F_B = x[1]  # molar flow rate of H2
    F_C = x[2]  # molar flow rate of Metanol
    F_D = x[3]  # molar flow rate of Agua
    F_E = x[4]  # molar flow rate of CO
    F_F = x[5]  # molar flow rate of CH4
    # compute total molar flow rate
    F_tot = F_A + F_B + F_C + F_D + F_E + F_F
    # compute volumetric flow rate in m^3/s
    R = Reactor.R  # Bar*L/mol*K
    T = Reactor.T_0  # isothermal operation
    Z = 1  # Acentric Factor
    vflow = (F_tot / Reactor.F_tot0) * (Reactor.P_0 / Reactor.P_0) * (
            Reactor.T_0 / T) * Reactor.vflow_0  # volumetric flow rate
    # compute partial pressures in gas phase
    P_A = (F_A * Z * R * T) / vflow  # bar
    P_B = (F_B * Z * R * T) / vflow  # bar
    P_C = (F_C * Z * R * T) / vflow  # bar
    P_D = (F_D * Z * R * T) / vflow  # bar
    P_E = (F_E * Z * R * T) / vflow  # bar
    P_F = (F_F * Z * R * T) / vflow  # bar
    # Compute reaction rate
    # Inhibition Parameter
    inh = (1 + (Rate.K_ads_co2 * P_A) + np.sqrt(Rate.K_ads_h2 * P_B)) ** 2
    # Rate law 3
    num_r3 = (np.sqrt(P_A * P_B) * (
            1 - ((P_F * (P_D ** 2)) / (P_A * (P_B ** 4) * Rate.k_eq_3))))  # Numerador velocidad de reacción 3
    r3F = (Rate.k_vel_3 * num_r3) / inh  # mol/kg cat*s
    return r3F


con1 = {'type': 'ineq', 'fun': constraint_rC}
con2 = {'type': 'ineq', 'fun': constraint_rCO2}
con3 = {'type': 'ineq', 'fun': constraint_rCH4}
b = (0, np.inf)
bnds = (b, b, b, b, b, b)
cons = [con1, con2, con3]
x0 = [1.10, 17.05, 0.542, 2.54, 1.705, 0.144]  # Values from PFR
sol = minimize(min_error_Balance, x0, args=(Reactor, Rate), method='SLSQP', bounds=bnds, constraints=cons, tol=1e-12)

# Selectividad Global
S_1 = np.divide(sol.x[2], np.add(sol.x[4],sol.x[5]))
#Rendimiento 1 Global
Y_1_1 = np.divide(sol.x[2],Reactor.F_A0-sol.x[0])

Sel_metanol = sol.x[2] / (Reactor.F_A0 - sol.x[0])
Conv_CO2 = (Reactor.F_A0 - sol.x[0]) / Reactor.F_A0
Prod_Metanol = (sol.x[2]) * (32.04 / 1000)

print("Selectividad Metanol libro rCH3OH/rCH4+rCO", S_1)
print("Rendimiento libro", Y_1_1)
print("Selectividad Metanol paper", Sel_metanol)
print("Conversion CO2 paper", Conv_CO2)
print("Produccion Metanol por lote", Prod_Metanol, "Kg/s")

stop = timeit.default_timer()
print('Time: ', stop - start)
