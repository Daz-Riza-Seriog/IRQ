# Code made for Sergio Andrés Díaz Ariza
# 10 September 2022
# License MIT
# IRQ: Python Program-Assigment 1


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from Parametros import Parametros
import timeit

start = timeit.default_timer()
sns.set()

Rate = Parametros.set_Rate()  # llamamos las K'as de las Reaciones
Reactor = Parametros.set_Reactor()  # llamamos los Parametros del reactor
ini_state = Parametros.set_initial_state_flux()  # llamamos los Estados iniciales del reactor
Est = Parametros.Estequiometria()  # llamamos Realaciones de Estequiometria


# Definimos la Ecuacion de Volumen en terminos de X
def Vol_CSTR(x):
    num_ = (Reactor.vflow_0 * (1 + Est.Eps * x) * x)
    den_ = (Rate.k_s * (1 - x))
    V = num_/den_
    return V


# Solucionamos la Ecuacion de Volumen asignandole una conversion
V_CSTR_09 = Vol_CSTR(Reactor.X)
print("Volumen para CSTR a X=", Reactor.X, ":", V_CSTR_09, "[m^3]")
print("Damkohler a X=", Reactor.X, ":", (V_CSTR_09*Rate.k_s/Reactor.vflow_0), "[m^3]")

conv = np.linspace(0, Reactor.X, 1000)
V_Cstr = [Vol_CSTR(x) for x in conv]
Da = [((x*Rate.k_s)/Reactor.vflow_0) for x in V_Cstr]


class plot_results:
    plt.figure(1)
    plt.title("CSTR\nConversion vs Volume")
    plt.plot(conv, V_Cstr)
    plt.xlabel('X', labelpad=15, fontsize=13)
    plt.ylabel('$[m{^3}]$', labelpad=8, fontsize=12)
    plt.tight_layout()

    plt.figure(2)
    plt.title("CSTR\nDamkohler vs Volume")
    plt.plot(Da, conv)
    plt.xlabel('$Da$', labelpad=15, fontsize=13)
    plt.ylabel('$X$', labelpad=8, fontsize=12)
    plt.tight_layout()

    plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)
