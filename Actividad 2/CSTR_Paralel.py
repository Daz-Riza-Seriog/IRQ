# Code made for Sergio Andrés Díaz Ariza
# 10 September 2022
# License MIT
# IRQ: Python Program-Assigment 1


import seaborn as sns
from Parametros import Parametros
from CSTR import V_CSTR_09
import timeit

start = timeit.default_timer()
sns.set()

Rate = Parametros.set_Rate()  # llamamos las K'as de las Reaciones
Reactor = Parametros.set_Reactor()  # llamamos los Parametros del reactor
Est = Parametros.Estequiometria()  # llamamos Realaciones de Estequiometria

# Preguntamos al usuario cuantos reactores en serie quiere
# n_div = int(input('Enter number of CSTR in series that you want:\n (condition 2 =< n)  '))
n_div = 3
# Dividimos equitativamente el Volumen del CSTR hallado anteriormente en las n divisiones
V_n = V_CSTR_09 / n_div


# Definimos la Ecuacion de Volumen en terminos de X para Cambios en el Flujo de entrada y Volumenes iguales
def Vol_CSTR_flux(x, n, frac):
    num_ = (Reactor.vflow_0 * (1 + Est.Eps * x) * x) * (frac)
    den_ = (Rate.k_s * (1 - x))
    V = num_ / den_
    return V


# Solucionamos la Ecuacion de Volumen asignandole una conversion, numero de reactores y fracciones de flujo
V_par_05 = Vol_CSTR_flux(Reactor.X, n_div, 0.5)
V_par_02 = Vol_CSTR_flux(Reactor.X, n_div, 0.2)
V_par_01 = Vol_CSTR_flux(Reactor.X, n_div, 0.3)

print("Volumen para CSTR paralelo 0.5F_A0:", V_par_05)
print("Volumen para CSTR paralelo 0.2F_A0:", V_par_02)
print("Volumen para CSTR paralelo 0.1F_A0:", V_par_01)


# Definimos la Ecuacion de Volumen en terminos de X para Cambios en los Volumenes de Recotres y FA0 igual
def Vol_CSTR_V(x, n, frac):
    num_ = (Reactor.vflow_0 * (1 + Est.Eps * x) * x)
    den_ = (Rate.k_s * (1 - x)) * (n * frac)
    V = num_ / den_
    return V


# Solucionamos la Ecuacion de Volumen asignandole una conversion, numero de reactores y fracciones de flujo
V_par_05_v = Vol_CSTR_V(Reactor.X, n_div, 0.5)
V_par_02_v = Vol_CSTR_V(Reactor.X, n_div, 0.2)
V_par_01_v = Vol_CSTR_V(Reactor.X, n_div, 0.3)

print("Volumen para CSTR paralelo con supuesto 0.5V.total :", V_par_05_v)
print("Volumen para CSTR paralelo supuesto 0.2V.total:", V_par_02_v)
print("Volumen para CSTR paralelo supuesto 0.3V.total:", V_par_01_v)


# Definimos la Ecuacion de Volumen en terminos de X para Cambios en los Volumenes de Recotres y FA0 igual
def Vol_CSTR_FV(x, frac_v, frac_f):
    num_ = (Reactor.vflow_0 * (1 + Est.Eps * x) * x) * frac_f
    den_ = (Rate.k_s * (1 - x)) * frac_v
    V = num_ / den_
    return V


# Solucionamos la Ecuacion de Volumen asignandole una conversion, numero de reactores y fracciones de flujo
V_par_05_vf = Vol_CSTR_FV(Reactor.X, 0.5, 0.3)
V_par_02_vf = Vol_CSTR_FV(Reactor.X, 0.2, 0.2)
V_par_01_vf = Vol_CSTR_FV(Reactor.X, 0.3, 0.5)

print("Volumen para CSTR paralelo con supuesto 0.5V.total y 0.1 FA0:", V_par_05_vf)
print("Volumen para CSTR paralelo supuesto 0.2V.total y 0.2 FA0:", V_par_02_vf)
print("Volumen para CSTR paralelo supuesto 0.1V.total y 0.5 FA0:", V_par_01_vf)

stop = timeit.default_timer()
print('Time: ', stop - start)
