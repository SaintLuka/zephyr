import numpy as np
import matplotlib.pyplot as plt
import os

from csv import CsvFile

exact_vars = ["x", "rho_exact", "u_exact", "e_exact", "p_exact"]
all_vars = ["x", "rho_exact", "u_exact", "e_exact", "p_exact", "rho", "u", "e", "p"]


def limits(arr):
    arr_min = np.min(arr)
    arr_max = np.max(arr)
    spread = arr_max - arr_min

    margin = 0.15
    arr_min -= margin * spread
    arr_max += margin * spread
    return np.round(arr_min, 1), np.round(arr_max, 1)


def plot_graphic(data: CsvFile, filename: str, title: str = None, with_numeric: bool = True):
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

    axr.plot(X, R, color='black', label='exact')
    axu.plot(X, U, color='black', label='exact')
    axe.plot(X, E, color='black', label='exact')
    axp.plot(X, P, color='black', label='exact')

    if with_numeric:
        axr.plot(X, data.rho, label='Godunov', linewidth=1.2)
        axu.plot(X, data.u, label='Godunov', linewidth=1.2)
        axe.plot(X, data.e, label='Godunov', linewidth=1.2)
        axp.plot(X, data.p, label='Godunov', linewidth=1.2)

    if title is not None:
        plt.title(title)

    plt.legend()
    fig.tight_layout()

    plt.savefig(filename)
    plt.close(fig)


def plot_graphic_compare(data1: CsvFile, data2: CsvFile, filename: str):
    X = data1.x
    R = data1.rho_exact
    U = data1.u_exact
    E = data1.e_exact
    P = data1.p_exact

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

    axr.plot(X, R, color='black', label='exact')
    axu.plot(X, U, color='black', label='exact')
    axe.plot(X, E, color='black', label='exact')
    axp.plot(X, P, color='black', label='exact')

    axr.plot(X, data1.rho, label='Godunov', linewidth=1.2)
    axu.plot(X, data1.u, label='Godunov', linewidth=1.2)
    axe.plot(X, data1.e, label='Godunov', linewidth=1.2)
    axp.plot(X, data1.p, label='Godunov', linewidth=1.2)

    axr.plot(X, data2.rho, label='Godunov 2 order', linewidth=1.2)
    axu.plot(X, data2.u, label='Godunov 2 order', linewidth=1.2)
    axe.plot(X, data2.e, label='Godunov 2 order', linewidth=1.2)
    axp.plot(X, data2.p, label='Godunov 2 order', linewidth=1.2)

    plt.legend()
    fig.tight_layout()

    plt.savefig(filename)
    plt.close(fig)


dir = "output"
# plot_graphic_compare(CsvFile('output/toro_test2_1_order.csv', all_vars),
#                      CsvFile('output/toro_test2_2_order.csv', all_vars),
#                      'graphics/toro_test2_compare')

plot_graphic_compare(CsvFile('output/sod_test_1_order.csv', all_vars),
                     CsvFile('output/sod_test_2_order.csv', all_vars),
                     'graphics/sod_test_compare')

# for filename in os.listdir(dir):
#     numeric = True
#     if filename.find("exact") != -1:
#         numeric = False
#
#     if numeric:
#         data = CsvFile(os.path.join(dir, filename), ["x", "rho_exact", "u_exact", "e_exact", "p_exact", "rho", "u", "e", "p"])
#     else:
#         data = CsvFile(os.path.join(dir, filename), ["x", "rho_exact", "u_exact", "e_exact", "p_exact"])
#
#     substr = filename[filename.find("toro_test") + len("toro_test") + 1:-4]
#     test_num, g1, g2 = substr.split('_')
#     png_name = os.path.join("graphics", filename[:-4].replace('.', ''))
#     try:
#         plot_graphic(data, png_name, with_numeric=numeric)
#     except Exception as e:
#         print(f'Error on filename: {filename}')
