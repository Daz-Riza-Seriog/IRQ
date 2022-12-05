# Code made by Sergio Andrés Díaz Ariza
# 10 September 2022
# License MIT
# IRQ: Python Program-Assigment 4


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit
import timeit

start = timeit.default_timer()
sns.set()


# ARRAYS


def func_data(x, alpha, beta, k):
    return ((C_A0 ** (1 - alpha)) - ((3 * C_MHF[-1]) ** beta) * (1 - alpha) * k * x) ** (1 / (1 - alpha))


t = [7.856861177, 16.61115681, 29.74513538, 45.52882719, 59.56733857, 90.32246797, 120.202472, 150.545897, 180.0131827]
C_A0 = 0.093778802
C_MHF = [0.010138249, 0.020529954, 0.037511521, 0.049930876, 0.06235023, 0.070207373, 0.080345622, 0.08718894,
         0.092764977]
C_A_exp = [0.093778802, 0.081866359, 0.064631336, 0.049677419, 0.038525346, 0.025599078, 0.013940092, 0.008110599,
           0.003294931]

xdata = t
a, b, k = 1.1, 1.2, 0.121  # suggest values for iteration
ydata = [func_data(x, a, b, k) for x in xdata]
# popt->Optimal values for the parameters so that the sum of the squared residuals of f(xdata, *popt) - ydata is minimized.
# pcov->The estimated covariance of popt. The diagonals provide the variance of the parameter estimate.
popt, pcov = curve_fit(func_data, xdata, ydata, bounds=(0, [3.0, 3.0, 1.0]), p0=[1.3, 1.9, 0.02], maxfev=5000)
# Note that p0 -> remarque the zone of the solution, is need have has a good suggest
print(popt)

y_fix = [func_data(x, *popt) for x in xdata]
print(ydata)
print(y_fix)
plt.figure(1)
plt.plot(xdata, ydata, 'b-', label='data')
plt.plot(xdata, y_fix, 'g--', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
plt.show()

stop = timeit.default_timer()
print('Time: ', stop - start)
