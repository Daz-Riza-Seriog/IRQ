# Code made for Sergio Andrés Díaz Ariza
# 10 September 2022
# License MIT
# IRQ: Python Program-Assigment 1


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from PIL._imaging import display

from Parametros import Parametros
from scipy.optimize import fsolve
from CSTR import V_CSTR_09
import timeit

start = timeit.default_timer()
sns.set()

Rate = Parametros.set_Rate()  # llamamos las K'as de las Reaciones
Reactor = Parametros.set_Reactor()  # llamamos los Parametros del reactor
Est = Parametros.Estequiometria()  # llamamos Realaciones de Estequiometria

# Preguntamos al usuario cuantos divisiones al Volumen inicial de CSTR quiere hacer
n_div = 10
# Dividimos equitativamente el Volumen del CSTR hallado anteriormente en las n divisiones
# TODO modificar 10 por V_CSTR_09
V_n = V_CSTR_09 / n_div



def V_(x, x_1):
    num_v = ((x[0] - x_1) * (1 + Est.Eps * x[0]) * Reactor.vflow_0)
    den_v = (Rate.k_s * (1 - x[0]))
    #0= num_v / den_v -V
    # [((((x[0] - x_1) * (1 + Est.Eps * x[0]) * Reactor.vflow_0) / (Rate.k_s * (1 - x[0]))) - V_n)]
    return [(num_v / den_v) -V_n]

# Lista de Reactores
n_list = [x + 1 for x in np.arange(n_div)]

# Loop de la primera conversion hacia la ultima
indx = 0
indx2 = 0
conv_ = 0  # Converesion inicial a la entrada del ciclo
x_list2 = np.zeros(n_div)
for x_n_1 in x_list2:

    indx2 = indx2
    arg = x_list2[indx2]

    X_n = fsolve(V_, 0.999, args=arg)
    if indx <= len(x_list2):
        x_list2[indx] = X_n
    else:
        break
    indx += 1
    indx2 = indx-1

Vol = [x*V_n for x in n_list]

d = {"Reactor": n_list, "Conversion": x_list2, "Volumen":Vol}
df = pd.DataFrame(d)
print(df)



class plot_results:
    plt.figure(1)
    plt.title(f"CSTR-series\nConversion vs ${n_div}$-Reactores")
    plt.plot(n_list, x_list2, marker="o", markersize=10, markeredgecolor="red", markerfacecolor="green")
    plt.xticks(np.arange(1, n_div, 1))
    plt.xlabel('$n$-Reactores', labelpad=15, fontsize=13)
    plt.ylabel('$X$', labelpad=8, fontsize=12)
    plt.tight_layout()

    plt.figure(2)
    plt.title(f"CSTR-series\nConversion vs Volumen")
    plt.plot(x_list2, Vol ,marker="o", markersize=10, markeredgecolor="red", markerfacecolor="green")
    plt.xlabel('$X$', labelpad=15, fontsize=13)
    plt.ylabel('Volumen $[m^3]$', labelpad=8, fontsize=12)
    plt.tight_layout()

    plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)
