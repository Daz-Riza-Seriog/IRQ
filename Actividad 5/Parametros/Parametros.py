# Code made by Sergio Andrés Díaz Ariza
# 04 November 2022
# License MIT
# IRQ: Python Program-Assigment 5

import numpy as np


class set_Reactor:
    Wf = 10  # kg de catalizador, este es un supuesto
    V_cstr = 100  # m3 supuesto volumen reactor

    # Condiciones a la entrada del reactor
    P0 = 101325  # [Pa]
    T0 = 673  # [K]
    Ta0 = 723  # [K] Temperatura inicial refrigerante
    # Constante universal de los gases
    Rg = 8.314  # Pa m3 /mol K

    # FLUJOS MOLARES
    F_cool = 2  # [mol/s] Flujo del refrigerante

    FA0 = 316.05  # Flujo de alimentación inicial de CO2 [mol/h] para una producción de 74.07m3/dia
    FA0 = FA0 / 3600  # Flujo de alimentación inicial de CO2 [mol/s]
    FB0 = 4 * FA0  # Flujo de alimentación incial de Hidrogeno [mol/s]
    FC0 = 0  # Flujo de alimentación incial de Metano [mol/s]
    FD0 = 0  # Flujo de alimentación incial de Agua [mol/s]
    FE0 = 11*FA0
    FT0 = FA0 + FB0 + FC0 + FD0 + FE0  # Flujo total de entrada [mol/s]

    # Molar Fraction in the inlet of CO2
    y_co2_o = FA0 / FT0
    # Molar Fraction in the inlet of H2
    y_h2_o = 1 - y_co2_o
    # partial pressure of A in inlet stream
    p_A0 = y_co2_o * P0  # Bar
    # partial pressure of I in inlet stream
    p_B0 = y_h2_o * P0  # Bar
    # Concentrations
    C_A0 = p_A0 / (Rg * T0)
    C_B0 = p_B0 / (Rg * T0)
    # Voumetric Initial Flow
    vflow_0 = FA0 / C_A0  # L/s

    # Pesos Moleculares
    PM_A = 44.01  # [g/mol] CO2
    PM_B = 2 * 1.00784  # [g/mol] H2
    PM_C = 16.04  # [g/mol] CH4
    PM_D = 18.01529  # [g/mol] H2O
    PM_E = 28.01  # [g/mol] CO


    # Parámetros caída de presión
    D = 0.32  # Diámetro reactor [m]
    Ac = np.pi * (D / 2) ** 2  # Área transversal [m^2]
    mu = 2.64e-5  # Visocisidad de la mezlca Pa s [kg/m*s]
    ro = 0.1756  # Densidad inicial [kg/m3];
    Dp = 0.005  # Diametro de partícula [m]
    f = 0.4  # Porosidad
    rc = 760  # kg/m3 densidad del catalizador
    Po = 101325  # Pa

    # G: superficial mass
    # velocity(kg / m2 / s)
    # % f: porosity
    # Ac: Crosssectional area(m2)
    # rc: catalyst density(kg / m3)
    # Po: entrance pressure(Pa)
    # ro: gas density(kg / m3)
    # Dp: particle diamenter(m)
    # mu: gas viscosity(kg / m / s)

Reactor = set_Reactor()


class set_Rate:
   # set reaction rate law parameters
    K_eq = 1.470e+02  # K de equilibrio [1/bar^2]

    # Párametros necesarios para la solución de las velocidad de reacción.

    # Calores de Absorcion dados en Articulo
    QAb_H2 = 52.0  # Calor H2[kJ/mol]
    QAb_CO2 = 9.72  # Calor CO2[kJ/mol]
    QAb_H2O = 14.5  # Caloe H2O[kJ/mol]
    QAb_CO = 40.6  # Caloe H2O[kJ/mol]

    # Factores frecuencia y Energias de Activacion para cada Rxn
    kf_1 = 1.14e8  # Constante de frecuencia CO2 methanation[mol min-1 g-1]
    Ea_1 = 110  # Energía de activación CO2 methanation[kJ/mol]
    kf_2 = 1.78e6  # Constante de frecuencia RWGS[mol min-1 g-1]
    Ea_2 = 97.1  # Energía de activación RWGS[kJ/mol]
    kf_3 = 2.23e8  # Constante de frecuencia CO methanation[mol min-1 g-1]
    Ea_3 = 97.3  # Energía de activación CO methanation[kJ/mol]

    R = 0.008314  # Constante de los gases [kJ/K *mol]

    KAb_H2 = 5.20E-5  # bar^-1 Constante de obtenida del artículo
    KAb_CO2 = 1.07  # bar^-1
    KAb_H2O = 6.09E-1  # bar^-1
    KAb_CO = 2.39E-3  # bar^-1

    # # CONSTANTE DE EQUILIBRIO ESTANDAR PARA CADA REACCION #
    T_est = 298.15  # Temperatura estandar [K]
    G_est_1 = -113.5  # Delta G estantar de reacción 1 [KJ/mol]
    G_est_2 = 28.62  # Delta G estantar de reacción 2 [KJ/mol]
    G_est_3 = -142.12  # Delta G estantar de reacción 3 [KJ/mol]

    k_eq_est1 = np.exp(-G_est_1 / (R * T_est))
    k_eq_est2 = np.exp(-G_est_2 / (R * T_est))
    k_eq_est3 = np.exp(-G_est_3 / (R * T_est))

    # ENTALPIAS PARA CADA REACCION
    H_rxn_1_est = -165  # Entalpia de reacción 1 en [KJ/mol]
    H_rxn_2_est = 41  # Entalpia de reacción 2 en [KJ/mol]
    H_rxn_3_est = -206  # Entalpia de reacción 3 en [KJ/mol]

    # # Constantes Equilibrio
    # k_eq_1 = Rate.k_eq_est1 * np.exp((Rate.H_rxn_1_est / R) * ((1 / Rate.T_est) - (1 / Reactor.T0)))  # Constante de equilibrio 1
    # k_eq_2 = Rate.k_eq_est2 * np.exp((Rate.H_rxn_2_est / R) * ((1 / Rate.T_est) - (1 / Reactor.T0)))  # Constante de equilibrio 2
    # k_eq_3 = Rate.k_eq_est3 * np.exp((Rate.H_rxn_3_est / R) * ((1 / Rate.T_est) - (1 / Reactor.T0)))  # Constante de equilibrio 3

    # entalpías a cualquier temperatura de cada reacción. Parte que no varía.
    # CPs para cada comoponente como promedio de 623-723 K
    Cp_CO2 = 48.7 / 1000  # [Kj/mol*K] Aspen
    Cp_CO = 30.95 / 1000  # [Kj/mol*K] Aspen
    Cp_H2 = 29.363 / 1000  # [Kj/mol*K] Aspen
    Cp_H2O = 37.17 / 1000  # [Kj/mol*K] Aspen
    Cp_CH4 = 56.6 / 1000  # [Kj/mol*K] Aspen

    temp_ref = 25 + 273.15
    D_Cp_1 = (1 / 1) * Cp_CH4 + (2 / 1) * Cp_H2O - (4 / 1) * Cp_H2 - Cp_CO2
    D_Cp_2 = (1 / 1) * Cp_CO + (1 / 1) * Cp_H2O - (1 / 1) * Cp_H2 - Cp_CO2
    D_Cp_3 = (1 / 1) * Cp_CH4 + (1 / 1) * Cp_H2O - (3 / 1) * Cp_H2 - Cp_CO

    # Global Coefficient of Transfer Heat -TODO toca hallarlo esto es inventado
    Ua = 8  # [kJ/s*m^3*K]
    # Cp para refrigerante
    Cp_cool = 107.5/1000  # [Kj/mol*K] nitrato de sodio

Rate = set_Rate()

class set_initial_state_flux_drop:
    __metaclass__ = Rate
    __metaclass__ = Reactor
    x_0 = np.empty(11)  # allocate memory
    x_0[0] = Reactor.FA0  # CO2
    x_0[1] = Reactor.FB0  # H2
    x_0[2] = np.finfo(float).eps  # CH4
    x_0[3] = np.finfo(float).eps  # H2O
    x_0[4] = Reactor.FE0  # CO
    x_0[5] = Reactor.P0  # pressure
    x_0[6] = Reactor.Ta0  # Initial Temperature of Fluid Cooler
    x_0[7] = Reactor.T0  #Initial Temperature of Reactor
    x_0[8] = np.finfo(float).eps # r1
    x_0[9] = np.finfo(float).eps # r2
    x_0[10] = np.finfo(float).eps # r3


class set_initial_state_flux_drop_adia:
    __metaclass__ = Rate
    __metaclass__ = Reactor
    x_0 = np.empty(10)  # allocate memory
    x_0[0] = Reactor.FA0  # molar flow rate of A
    x_0[1] = Reactor.FB0  # molar flow rate of B
    x_0[2] = np.finfo(float).eps  # molar flow rate of C
    x_0[3] = np.finfo(float).eps  # (Rate.k_eq_3 * (np.sqrt(Reactor.p_A0 * Reactor.p_B0))) / Rate.inh_0
    x_0[4] = Reactor.FE0  # (Rate.k_eq_3 * (np.sqrt(Reactor.p_A0 * Reactor.p_B0))) / Rate.inh_0
    x_0[5] = Reactor.P0  # pressure
    x_0[6] = Reactor.T0  # Initial Temperature of Reactor
    x_0[7] = np.finfo(float).eps
    x_0[8] = np.finfo(float).eps
    x_0[9] = np.finfo(float).eps

class set_initial_state_flux_cstr:
    __metaclass__ = Rate
    __metaclass__ = Reactor
    x_0 = np.empty(5)  # allocate memory
    x_0[0] = Reactor.FA0  # molar flow rate of A
    x_0[1] = Reactor.FB0  # molar flow rate of B
    x_0[2] = np.finfo(float).eps  # molar flow rate of C
    x_0[3] = np.finfo(float).eps  # (Rate.k_eq_3 * (np.sqrt(Reactor.p_A0 * Reactor.p_B0))) / Rate.inh_0
    x_0[4] = Reactor.FE0  # (Rate.k_eq_3 * (np.sqrt(Reactor.p_A0 * Reactor.p_B0))) / Rate.inh_0