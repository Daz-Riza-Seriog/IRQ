# Code made for Sergio Andrés Díaz Ariza
# 17 September 2022
# License MIT
# IRQ: Python Program-Assigment 1

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import timeit

# import the classes for the functions
from Arrhenius import Arrhenius
from Velocity_Rxn import Velocity_Rxn

start = timeit.default_timer()
sns.set()

## Parameters ##
# Initial Concentrations
C_H2SO4_o = 1  # [M] Initial Concentration
C_HMF_o = 0.1  # [M] Initial Concentration
C_H2O_o = 0.6  # [M]
C_FA_o = 0  # [M]
C_LA_o = 0  # [M]

# Activation Energies
E_1 = 110.5  # [KJ/mol]
E_2 = 111  # [KJ/mol]
R = 0.008314  # [KJ/mol*K] Constant ideal gas
T_r = 140 + 273.15  # [K] Temperature of reference in Arrhenius expression
X_ = 0.95  # Conversion percentage

# Parameters of Molar Balance
Theta_H2SO4 = C_H2SO4_o / C_HMF_o  # [mol H2SO4/mol HMF]
Theta_H2O = C_H2O_o / C_HMF_o  # [mol H2O/mol HMF]
Theta_FA = C_FA_o / C_HMF_o  # [mol FA/mol HMF]
Theta_LA = C_LA_o / C_HMF_o  # [mol LA/mol HMF]

F_HMF_o = 15  # [mol/s] Initial Concentration obtain from demand of LA

# Initial Flux molar
F_H2SO4_o = Theta_H2SO4 * F_HMF_o  # [M] Initial Concentration
F_H2O_o = Theta_H2O * F_HMF_o  # [M]
F_FA_o = Theta_FA * F_HMF_o  # [M]
F_LA_o = Theta_LA * F_HMF_o  # [M]

Arrh = Arrhenius.Arrhenius()
Vel_r = Velocity_Rxn.Velocity_Rxn()

# lists of Temperatures
T_ = np.arange(98, 182, 2)  # [C]
T = T_ + 273.15  # [K]
T_inv = 1 / T_

# Evaluate in Arrhenius for obtain K_1 for some T
K_1 = Arrh.K1_H(T, T_r, E_1, R)
r_1_ = Vel_r.r_1(K_1, C_HMF_o, Theta_H2SO4, X_)  # [M/s]
r_1 = r_1_ * (1000 / 60)  # [mol/m^3*s]
ln_r_1 = np.log(r_1)

# Crate a DataFrame of the data
df = pd.DataFrame(list(zip(T_, T, T_inv, K_1, r_1_, r_1, ln_r_1)),
                  columns=['T [C]', 'T [K]', '1/T [1/K]', 'K_1', 'r_1 [M/min]', 'r_1 [mol/m^3*s]', 'ln(r_1)'])

# Plot the Graphics r_MHF vs T
plt.figure(1)
plt.plot(T_, r_1)
plt.title("Velocity Law vs Temperature", fontsize=16)
plt.ylabel("$r_{HMF}$ $[mol/m^3*s]$", fontsize=14)
plt.xlabel("T $[\circ C]$ ", fontsize=14)

plt.figure(2)
plt.plot(T_inv, ln_r_1)
plt.title("Natural Logarithm Velocity Law vs Temperature", fontsize=16)
plt.ylabel("$\ln{r_{HMF}}$ $[adim]$", fontsize=14)
plt.xlabel("1/T $[\circ C^{-1}]$ ", fontsize=14)
df.plot.scatter('T [C]', 'r_1 [mol/m^3*s]')

plt.show()

stop = timeit.default_timer()
print('Time: ', stop - start)
