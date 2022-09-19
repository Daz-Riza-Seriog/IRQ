# Code made for Sergio Andrés Díaz Ariza
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
    C_A = x[0]  # molar flow rate of isopropanol
    C_B = x[1]  # molar flow rate of acetona
    C_C = x[2]  # molar flow rate of hidrogeno
    X = x[3]  # molar flow rate of isopropanol

    C_A = Est.C_A0 * (1 - X)
    C_B = Est.C_A0 * (Est.Theta_B0 + Est.b_a * X)
    C_C = Est.C_A0 * (Est.Theta_C0 + Est.c_a * X)
    C_I = Est.C_A0 * Est.Theta_I0

    # compute reaction rate
    r_R = Rate.k_s * C_A

    # Next, evaluate derivative functions.

    # mole balance on A
    dCA_dt = -r_R
    # mole balance on B
    dCB_dt = r_R
    # mole balance on C
    dCC_dt = r_R
    # conversion balance
    dXA_dt = r_R/Est.C_A0

    f = [dCA_dt,dCB_dt,dCC_dt,dXA_dt]

    return f


# Integramos (1+epsX)/(1-X)
def Int_Mole(x):
    f = 1 / (1 - x)
    return f


# Solucionamos para una lista de Conversiones de 0-0.9
conv = np.linspace(0, Reactor.X, 1000)
sol_int = [integrate.quad(Int_Mole, 0, x) for x in conv]
# Hallamos tiempo de residencia en Batch
t_batch = [i[0] * (1 / Rate.k_s) for i in sol_int]
tim_batch_09 = t_batch[-1]  # Tiempo para una Conversion de X=0.9
V_batch_09 = (Reactor.F_A0 / Reactor.C_A0) * tim_batch_09  # Volumen para Batch X=0.9
V_batch_09_25 = (Reactor.F_A0 / Reactor.C_A0) * (tim_batch_09+(2.5*3600))  # Volumen para Batch X=0.9 considerando Tiempos de calentamiento,vaciado y llenado

# Hallamos lista de Volumenes para una Conversion dada
V_batch = [((Reactor.F_A0 / Reactor.C_A0) * x) for x in t_batch]
Da_ = [(Rate.k_s*x) for x in t_batch]


print("Volumen para BATCH a X=", Reactor.X, ":", V_batch_09, "[m^3]")
print("Volumen para BATCH a X=", Reactor.X, ", con t= 2.5 horas de operacion:", V_batch_09_25, "[m^3]")
print("Damkohler para BATCH Da:", (Rate.k_s * tim_batch_09))



t_eval = np.linspace(0, Reactor.t_batch, 1000)
sol = solve_ivp(BATCH_calc, [0, Reactor.t_batch], ini_state.x_0, args=(Rate, Reactor, Est), t_eval=t_eval)

r_A = [(Rate.k_s*i) for i in sol.y[0]]

class plot_results:
    plt.figure(1)
    plt.title("BATCH\nTiempo vs Conversion")
    plt.plot(sol.t,sol.y[3])
    plt.xlabel('$[s]$', labelpad=15, fontsize=13)
    plt.ylabel('X', labelpad=8, fontsize=12)
    plt.tight_layout()

    plt.figure(2)
    plt.title("BATCH\nTime vs Conversion")
    plt.plot(sol.t, sol.y[0], label="isopropanol")
    plt.plot(sol.t, sol.y[1], label="acetona")
    plt.plot(sol.t, sol.y[2], label="hidrogeno")
    plt.xlabel('time $[s]$', labelpad=15, fontsize=13)
    plt.legend()
    plt.ylabel('Concentracion $C_i$ $[mol/m^3]$', labelpad=8, fontsize=12)
    plt.tight_layout()

    plt.figure(3)
    plt.title("BATCH\n Time vs $r_A$")
    plt.plot(sol.t, r_A)
    plt.xlabel('time $[s]$', labelpad=15, fontsize=13)
    plt.legend()
    plt.ylabel('$r_A$ $[mol/m^3*s]$', labelpad=8, fontsize=12)
    plt.tight_layout()

    plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)
