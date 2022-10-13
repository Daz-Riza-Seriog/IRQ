# Code made by Sergio Andrés Díaz Ariza
# 30 September 2022
# License MIT
# IRQ: Python Program-Assigment 1

import numpy as np


class set_Reactor:
    # Supose a Convertion
    X = 0.9
    # Supose a Initial Volume of CSTR
    V_cstr = 4  # m^3
    # Supose a Initial Volume of Reactor of PFR- Uncomment whe use it
    V_tot_su_pfr = 10  # Kg
    # Supose a Initial Volume of Reactor of PFR- Uncomment whe use it
    V_tot_su = 5  # Kg
    # Supose a Initial Time for Reactor Batch
    t_batch = 100  # [s]
    # 2.750

    ## VALORES INICIALES DADOS POR EL PAPER ##
    # Constant Ideal Gas for work in stoichiometry
    R = 8.3144626181532e-5  # Bar*m^3/K*mol
    # total pressure of inlet stream
    P_0 = 40  # in Bar initial
    # inlet temperature
    T_0 = 300 + 273.15  # in K
    # inlet volumetric flow rate of CO2
    F_A0 = 3.4936002609389596  # Valor dado en el trabajo anterior derivado de la meta a producir mol/s
    # molar flow rate of H2 in feed stream
    F_B0 = 6 * F_A0  # mol/s
    # molar flow rate of Metanol in feed stream
    F_C0 = 0  # mol/s
    # molar flow rate of Agua in feed stream
    F_D0 = 0  # mol/s
    # molar flow rate of CO in feed stream
    F_E0 = 0  # mol/s
    # molar flow rate of CH4 in feed stream
    F_F0 = 0  # mol/s
    # total molar flow rate of feed stream in mol/s
    F_tot0 = F_A0 + F_B0 + F_C0 + F_D0 + F_E0 + F_F0
    # Molar Fraction in the inlet of CO2
    y_co2_o = F_A0 / F_tot0

    # Molar Fraction in the inlet of H2
    y_h2_o = 1 - y_co2_o
    # partial pressure of A in inlet stream
    p_A0 = y_co2_o * P_0  # Bar
    # partial pressure of I in inlet stream
    p_B0 = y_h2_o * P_0  # Bar
    # Concentrations
    C_A0 = p_A0 / (R * T_0)
    C_B0 = p_B0 / (R * T_0)
    # Voumetric Initial Flow
    vflow_0 = F_A0 / C_A0  # L/s

    ## VALORES PARA EL LECHO-!!

    # AREA DEL REACTOR
    di = 0.70  # Diametro interno [m]
    area = np.pi * (di ** 2) / 4  # Area de los tubos #LE FALTO PONER PI TERRIBLE ERROR!!!

    # Calculo de Gamma ESTO ES DEL TRABAJO ANTERIOR CON APROXIMACIONES
    PM_A = 44.01 / 1000  # Peso molecular del co2 en [kg/mol]
    PM_B = 2.06 / 1000  # Peso molecular del h2 en [kg/mol]
    F_MASS = PM_A * F_A0 + PM_B * F_B0  # Flujo másico [kg/s]
    gamma = F_MASS / area  # kg/s*m2
    # Aproximación de la viscosidad del gas
    mu_h2 = 8.76e-7  # kg/m s
    mu_co2 = 1.48e-6  # kg/m s
    mu = (y_co2_o * mu_co2) + (y_h2_o * mu_h2)  # kg/m s
    # Aproximación de la densidad del gas
    PM_prom = (y_co2_o * PM_A) + (y_h2_o * PM_B)  # Peso molecular promedio de la mezcla [kg/mol]
    ro = (P_0 * PM_prom) / (T_0 * R)  # kg/m^3 Usando gas ideal
    # Diametro de particula
    Dp = 0.036  # m
    # Calculo de la porosidad
    f = 0.95  # Correlación para hallar la porosidad
    rc = 8780  # kg/m3 Densidad del catalizador #TODO-Aqui piden densidad particulas de catlizador-Es lo mimso?

    # Ergun constant beta_0 for model of pressure
    # drop across packed bed, assumed constant viscosity
    var1 = (150 * (1 - f) * mu) / Dp
    var2 = (gamma * (1 - f)) / (ro * Dp * (f ** 3))  # Here the Value from matlab is different var2=34.7669
    beta_0 = var2 * (var1 + 1.75 * gamma)

    # Coeficiente de transferencia de masa:
    K = 9.91496e-2  # [m3/kg cat*s]  # ms^-1


Reactor = set_Reactor()


class set_Rate:
    # set reaction rate law parameters
    R_1 = 8.314472e-3  # [KJ/mol*K] - We use that only for K´s determnation
    # TODO revisar las unidades con R, nuevamente los valores estan diferentes para concentraciones
    # CONSTANTE DE EQUILIBRIO ESTANDAR PARA CADA REACCION #
    T_est = 298.15  # Temperatura estandar [K]
    G_est_1 = -0.47  # Delta G estantar de reacción 1 [KJ/mol]
    G_est_2 = 23.4  # Delta G estantar de reacción 1 [KJ/mol]
    G_est_3 = -130.8  # Delta G estantar de reacción 1 [KJ/mol]
    k_eq_est1 = np.exp(-G_est_1 / (R_1 * T_est))
    k_eq_est2 = np.exp(-G_est_2 / (R_1 * T_est))
    k_eq_est3 = np.exp(-G_est_3 / (R_1 * T_est))

    # CONSTANTE DE EQUILIBRIO PARA CADA REACCION A TEMPERATURA DE OPERACION-CORRECCION VAN´T HOFF#
    H_rxn_1 = -49.479  # Entalpia de reacción 1 en [KJ/mol]
    H_rxn_2 = 41.161  # Entalpia de reacción 2 en [KJ/mol]
    H_rxn_3 = -164.367  # Entalpia de reacción 3 en [KJ/mol]
    k_eq_1 = k_eq_est1 * np.exp(
        (H_rxn_1 / R_1) * ((1 / T_est) - (1 / Reactor.T_0)))  # Constante de equilibrio 1
    k_eq_2 = k_eq_est2 * np.exp(
        (H_rxn_2 / R_1) * ((1 / T_est) - (1 / Reactor.T_0)))  # Constante de equilibrio 2
    k_eq_3 = k_eq_est3 * np.exp(
        (H_rxn_3 / R_1) * ((1 / T_est) - (1 / Reactor.T_0)))  # Constante de equilibrio 3

    # CONSTANTE DE ADSORCION CO2-H2 A TEMPERATURA DE OPERACION-VAN´T HOFF#
    K_ads_h2_300 = 0.76  # Constante de adsorción para el H2 [1/Bar]
    H_ads_h2_300 = -12.5  # [KJ/mol]
    K_ads_co2_300 = 0.79  # Constante de adsorción para el CO2 [1/Bar]
    H_ads_co2_300 = -25.9  # [KJ/mol]
    K_ads_h2 = K_ads_h2_300 * np.exp(
        (H_ads_h2_300 / R_1) * ((1 / Reactor.T_0) - (1 / Reactor.T_0)))  # Constante de adsorcion H2
    K_ads_co2 = K_ads_co2_300 * np.exp(
        (H_ads_co2_300 / R_1) * ((1 / Reactor.T_0) - (1 / Reactor.T_0)))  # Constante de adsorcion CO2

    # CONSTANTES DE VELOCIDAD PARA CADA REACCION REFERENCIA A 300 °C
    k_vel_1_ref = 6.9e-4  # Constante de velocidad de reacción 1 CH3OH [mol/s*bar^2k*kgcat]
    k_vel_2_ref = 1.8e-3  # Constante de velocidad de reacción 2 RWGS [mol/s*bar^1.5k*kgcat]
    k_vel_3_ref = 1.1e-4  # Constante de velocidad de reacción 3 Methanation [mol/s*bar^1k*kgcat]
    E_a_1 = 35.7  # [KJ/mol*K]
    E_a_2 = 54.5  # [KJ/mol*K]
    E_a_3 = 42.5  # [KJ/mol*K]

    # CONSTANTES DE VELOCIDAD PARA CADA REACCION REFERENCIA A T OPERACION
    k_vel_1 = k_vel_1_ref * np.exp((E_a_1 / R_1) * ((1 / T_est) - (1 / Reactor.T_0)))  # CH3OH [mol/s*bar^2k*kgcat]
    k_vel_2 = k_vel_2_ref * np.exp((E_a_2 / R_1) * ((1 / T_est) - (1 / Reactor.T_0)))  # RWGS [mol/s*bar^1.5k*kgcat]
    k_vel_3 = k_vel_3_ref * np.exp(
        (E_a_3 / R_1) * ((1 / T_est) - (1 / Reactor.T_0)))  # Methanation [mol/s*bar^1k*kgcat]

    # Inhibition Parameter Initial
    inh_0 = (1 + (K_ads_co2 * Reactor.p_A0) + np.sqrt(K_ads_h2 * Reactor.p_B0)) ** 2


Rate = set_Rate()


class Estequiometria():
    __metaclass__ = Reactor

    # Initial Concentrations
    C_A0 = Reactor.C_A0
    C_B0 = Reactor.C_B0
    C_C0 = Reactor.F_C0 / Reactor.vflow_0
    C_D0 = Reactor.F_D0 / Reactor.vflow_0
    C_E0 = Reactor.F_E0 / Reactor.vflow_0
    C_F0 = Reactor.F_F0 / Reactor.vflow_0

    # Thetas for each component
    Theta_B0 = C_B0 / C_A0  # Esta es para H2
    Theta_C0 = C_C0 / C_A0  # Esta es para Metanol
    Theta_D0 = C_D0 / C_A0  # Esta es para H20

    # Coeficientes Estequiometricos reaccion
    a = 1  # Coeficiente del CO2
    b = 3  # Coeficiente del H2
    c = 1  # Coeficiente del Metanol
    d = 1  # Coeficiente del H2O

    a_c = (a / c)  # Coeficiente estequiometrico a/c
    b_c = (b / c)  # Coeficiente estequiometrico b/a
    c_c = (c / c)  # Coeficiente estequiometrico c/a
    d_c = (d / c)  # Coeficiente estequiometrico d/a

    sigma = (a_c + b_c - d_c - 1)
    Eps = sigma * Reactor.y_co2_o


class set_initial_state_batch:
    __metaclass__ = Rate
    __metaclass__ = Reactor
    x_0 = np.empty(9)  # allocate memory
    x_0[0] = Reactor.C_A0  # molar concentration A
    x_0[1] = Reactor.C_B0  # molar concentration B
    x_0[2] = np.finfo(float).eps  # molar concentration C
    x_0[3] = np.finfo(float).eps  # molar concentration D
    x_0[4] = np.finfo(float).eps  # molar concentration E
    x_0[5] = np.finfo(float).eps  # molar concentration F
    x_0[6] = np.finfo(float).eps  # (Rate.k_eq_1 * (Reactor.p_A0 * (Reactor.p_B0 ** 2)) / Reactor.p_B0) / Rate.inh_0
    x_0[7] = np.finfo(float).eps  # (Rate.k_eq_2 * (Reactor.p_A0 * Reactor.p_B0) / np.sqrt(Reactor.p_B0)) / Rate.inh_0
    x_0[8] = np.finfo(float).eps  # (Rate.k_eq_3 * (np.sqrt(Reactor.p_A0 * Reactor.p_B0))) / Rate.inh_0
    # x_0[6] = 0  # conversion inicial


class set_initial_state_flux:
    __metaclass__ = Rate
    __metaclass__ = Reactor
    x_0 = np.empty(9)  # allocate memory
    x_0[0] = Reactor.F_A0  # molar flow rate of A
    x_0[1] = Reactor.F_B0  # molar flow rate of B
    x_0[2] = np.finfo(float).eps  # molar flow rate of C
    x_0[3] = np.finfo(float).eps  # molar flow rate of D
    x_0[4] = np.finfo(float).eps  # molar flow rate of E
    x_0[5] = np.finfo(float).eps  # molar flow rate of F
    x_0[6] = np.finfo(float).eps  # (Rate.k_eq_1 * (Reactor.p_A0 * (Reactor.p_B0 ** 2)) / Reactor.p_B0) / Rate.inh_0
    x_0[7] = np.finfo(float).eps  # (Rate.k_eq_2 * (Reactor.p_A0 * Reactor.p_B0) / np.sqrt(Reactor.p_B0)) / Rate.inh_0
    x_0[8] = np.finfo(float).eps  # (Rate.k_eq_3 * (np.sqrt(Reactor.p_A0 * Reactor.p_B0))) / Rate.inh_0


class set_initial_state_flux_drop:
    __metaclass__ = Rate
    __metaclass__ = Reactor
    x_0 = np.empty(10)  # allocate memory
    x_0[0] = Reactor.F_A0  # molar flow rate of A
    x_0[1] = Reactor.F_B0  # molar flow rate of B
    x_0[2] = np.finfo(float).eps  # molar flow rate of C
    x_0[3] = np.finfo(float).eps  # molar flow rate of D
    x_0[4] = np.finfo(float).eps  # molar flow rate of E
    x_0[5] = np.finfo(float).eps  # molar flow rate of F
    x_0[6] = np.finfo(float).eps  # (Rate.k_eq_1 * (Reactor.p_A0 * (Reactor.p_B0 ** 2)) / Reactor.p_B0) / Rate.inh_0
    x_0[7] = np.finfo(float).eps  # (Rate.k_eq_2 * (Reactor.p_A0 * Reactor.p_B0) / np.sqrt(Reactor.p_B0)) / Rate.inh_0
    x_0[8] = np.finfo(float).eps  # (Rate.k_eq_3 * (np.sqrt(Reactor.p_A0 * Reactor.p_B0))) / Rate.inh_0
    x_0[9] = Reactor.P_0  # pressure
