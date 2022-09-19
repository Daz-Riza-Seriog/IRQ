# Code made for Sergio Andrés Díaz Ariza
# 30 September 2022
# License MIT
# IRQ: Python Program-Assigment 1

import numpy as np


class set_Rate:
    # set reaction rate law parameters
    # forward surface reaction rate constant
    k_s = 129.81  # s^-1
    # k_s = 0.12981 # supuesta para ver como varian los CSTR en serie


class set_Reactor:
    # Supose a Convertion
    X = 0.9
    # Supose a Initial Volume of Reactor
    V_tot_su = 0.0043  # m^3
    # Supose a Initial Time for Reactor Batch
    t_batch = 0.03  # [s]

    # AREA DEL REACTOR
    l = 6.096  # Longitud de los tubos [m]
    tubos = 400  # Cantidad de tubos [m]
    di = 4.6 / 100  # Diametro interno [m]
    area = np.pi * (di ** 2) / 4  # Area de los tubos #LE FALTO PONER PI TERRIBLE ERROR!!!

    ## VALORES INICIALES DADOS POR EL PAPER ##
    # Constant Ideal Gas
    R = 8.31446261815324  # Pas*m^3/mol*K
    # total pressure of inlet stream
    P_0 = 3.3 * 1e5  # in Pas initial in 1.20 se modifico
    # inlet temperature
    T_0 = 723  # in K
    # Molar Fraction in the inlet of ISOPROPANOL
    y_iso_o = 0.7
    # partial pressure of A in inlet stream
    p_A0 = y_iso_o * P_0  # Pas
    # partial pressure of I in inlet stream
    p_I0 = (1 - y_iso_o) * P_0  # Pas
    # Concentrations
    C_A0 = p_A0 / (R * T_0)
    C_I0 = p_I0 / (R * T_0)
    # inlet volumetric flow rate
    F_A0 = 2.945075757575758 * (1000 / 3600)  # Valor dado en el trabajo anterior derivado de la meta a producir mol/s
    # F_A0 = 2945075.757575758  # Molar flux supose for CSTR in series changes
    vflow_0 = F_A0 / C_A0  # m^3/s
    # vflow_0 = 0.00723967
    # molar flow rate of B in feed stream
    F_B0 = 0  # mol/s
    # molar flow rate of B in feed stream
    F_C0 = 0  # mol/s
    # molar flow rate of I in feed stream
    F_I0 = C_I0 * vflow_0  # mol/s
    # total molar flow rate of feed stream in mol/s
    F_tot0 = F_A0 + F_B0 + F_C0 + F_I0

    ## VALORES PARA EL LECHO-!!SI ES QUE USAMOS!!
    # Densidad de catalizador
    d_cat = 2732  # kg/m^3 formado por cobre, zinc y cromo en un soporte de alúmina
    # diameter of catalyst particles
    D_p = 1.41 / 1000  # in m
    # void fraction of bed
    phi = 0.41  # Porosidad
    # cross sectional area
    A_c = 0.0016619025137490004  # in m ^ 2 "´POR DEFINIR"
    # density of solid catalyst phase
    rho_s = (1 - phi) * d_cat  # in lb cat/pie^3
    # density of gas phase at inlet conditions
    rho_0 = 0.12  # in Kg / m ^ 3
    # viscosity of gas pahse
    mu = 9.4e-6  # in Pa * s
    # superficial mass velocity
    gamma = (rho_0 * vflow_0) / A_c
    # Ergun constant beta_0 for model of pressure
    # drop across packed bed, assumed constant viscosity
    var1 = (150 * (1 - phi) * mu) / D_p
    var2 = (gamma * (1 - phi)) / (rho_0 * D_p * (phi ** 3))  # Here the Value from matlab is different var2=34.7669
    beta_0 = var2 * (var1 + 1.75 * gamma)
    # Coeficiente global de transferencia de masa
    K = 0.5  # ms^-1
    # Area sobre volumen
    a_difu = 4 / di


Rate = set_Rate()
Reactor = set_Reactor()


class Estequiometria():
    __metaclass__ = Reactor

    # Initial Concentrations
    C_A0 = Reactor.C_A0
    C_B0 = Reactor.F_B0 / Reactor.vflow_0
    C_C0 = Reactor.F_C0 / Reactor.vflow_0
    C_I0 = Reactor.C_I0

    # Thetas for each component
    Theta_B0 = C_B0 / C_A0
    Theta_C0 = C_C0 / C_A0
    Theta_I0 = C_I0 / C_A0

    # Coeficientes Estequiometricos
    a = 1
    b = 1
    c = 1

    b_a = (b / a)  # Coeficiente estequiometrico b/a
    c_a = (c / a)  # Coeficiente estequiometrico c/a

    sigma = (b_a + c_a - 1)
    Eps = sigma * Reactor.y_iso_o

class set_initial_state_batch:
    __metaclass__ = Rate
    __metaclass__ = Reactor
    x_0 = np.empty(4)  # allocate memory
    R = 8.314  # ideal gas constant in SI units
    x_0[0] = Reactor.C_A0
    x_0[1] = 0  # molar flow rate of B
    x_0[2] = 0  # molar flow rate of C
    x_0[3] = 0  # conversion inicial


class set_initial_state_flux:
    __metaclass__ = Rate
    __metaclass__ = Reactor
    x_0 = np.empty(3)  # allocate memory
    R = 8.314  # ideal gas constant in SI units
    x_0[0] = Reactor.F_A0
    x_0[1] = 0  # molar flow rate of B
    x_0[2] = 0  # molar flow rate of C


class set_initial_state_flux_drop:
    __metaclass__ = Rate
    __metaclass__ = Reactor
    x_0 = np.empty(4)  # allocate memory
    R = 8.314  # ideal gas constant in SI units
    x_0[0] = Reactor.F_A0
    x_0[1] = 0  # molar flow rate of B
    x_0[2] = 0  # molar flow rate of C
    x_0[3] = Reactor.P_0  # pressure
