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

dir = "ToroTest5"
names = ["Godunov", "Rusanov", "HLLC"]
data = []
for name in names:
    data.append(CsvFile(dir + "/results_" + name + ".csv", ["x", "rho_exact", "u_exact", "e_exact", "p_exact",
                                "rho", "u", "e", "p"]))

data = np.array(data)
X = data[0].x
R = data[0].rho_exact
U = data[0].u_exact
E = data[0].e_exact
P = data[0].p_exact

fig = plt.figure(dpi=400, figsize=(10.8, 6.4))

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

axr.plot(X, R, color='black', label='exact')
axu.plot(X, U, color='black', label='exact')
axe.plot(X, E, color='black', label='exact')
axp.plot(X, P, color='black', label='exact')

for arr, name in zip(data, names):
    axr.plot(X, arr.rho, label=name, linewidth=1.2)
    axu.plot(X, arr.u, label=name, linewidth=1.2)
    axe.plot(X, arr.e, label=name, linewidth=1.2)
    axp.plot(X, arr.p, label=name, linewidth=1.2)

fig.tight_layout()

plt.legend()

plt.savefig(dir + "/result.png")