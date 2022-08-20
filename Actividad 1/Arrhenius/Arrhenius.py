# Code made for Sergio Andrés Díaz Ariza
# 20 September 2022
# License MIT
# IRQ: Python Program-Assigment 1

import numpy as np

class Arrhenius:

    def K1_H(self, T, T_r, E_1, R):
        K1 = 0.340 * np.exp((-E_1 / R) * ((1 / T) - (1 / T_r)))
        return K1

    def K2_H(self, T, T_r, E_2, R):
        K2 = 0.117 * np.exp((-E_2 / R) * ((1 / T) - (1 / T_r)))
        return K2