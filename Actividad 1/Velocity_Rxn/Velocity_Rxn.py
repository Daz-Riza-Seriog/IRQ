# Code made for Sergio Andrés Díaz Ariza
# 20 September 2022
# License MIT
# IRQ: Python Program-Assigment 1


class Velocity_Rxn:

    def r_1(self, K1_H, C_HMF_o, Theta_H2SO4, X):
        r_1 = K1_H * ((C_HMF_o * (1 - X)) ** 0.88) * (Theta_H2SO4 * C_HMF_o) ** 1.38
        return r_1

    def r_2(self, K2_H, C_HMF_o, Theta_H2SO4, X):
        r_1 = K2_H * ((C_HMF_o * (1 - X)) ** 1.23) * (Theta_H2SO4 * C_HMF_o) ** 1.09
        return r_1

