"""
An example of reading .csv
"""
import numpy as np
import matplotlib.pyplot as plt

from csv import CsvFile


def limits(arr):
    arr_min = np.min(arr)
    arr_max = np.max(arr)
    spread = arr_max - arr_min

    margin = 0.15
    arr_min -= margin * spread
    arr_max += margin * spread
    return np.round(arr_min, 1), np.round(arr_max, 1)


data = CsvFile("results.csv", ["x", "rho_exact", "u_exact", "e_exact", "p_exact"])

X = data.x
R = data.rho_exact
U = data.u_exact
E = data.e_exact
P = data.p_exact


fig = plt.figure(dpi=140, figsize=(10.8, 6.4))

axr = fig.add_subplot(2, 2, 1)
axu = fig.add_subplot(2, 2, 2)
axe = fig.add_subplot(2, 2, 3)
axp = fig.add_subplot(2, 2, 4)

axr.set_title('Плотность')
axu.set_title('Скорость')
axe.set_title('Энергия')
axp.set_title('Давление')

x_min = np.round(np.min(X), 2)
x_max = np.round(np.max(X), 2)
axr.set_xlim(x_min, x_max)
axu.set_xlim(x_min, x_max)
axe.set_xlim(x_min, x_max)
axp.set_xlim(x_min, x_max)

axr.set_ylim(limits(R))
axu.set_ylim(limits(U))
axe.set_ylim(limits(E))
axp.set_ylim(limits(P))

axr.plot(X, R, color='black')
axu.plot(X, U, color='black')
axe.plot(X, E, color='black')
axp.plot(X, P, color='black')

fig.tight_layout()

plt.show()