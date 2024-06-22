"""
An example of reading .csv
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from csv import CsvFile




data = CsvFile("test2D.csv", ["x", "y", "rho"])
data.as2D()

X = data.x
Y = data.y
R = data.rho

fig = plt.figure(dpi=140, figsize=(5.2, 4.8))

ax = plt.gca()

ax.set_title('Плотность')

divider = make_axes_locatable(ax)
ax_cb = divider.append_axes("right", size="5%", pad=0.1)
fig.add_axes(ax_cb)

imshow = True

if imshow:
    options = {
        'cmap': 'jet',
        'origin':'lower',
        'extent' : [data.xmin, data.xmax,
                    data.ymin, data.ymax],
        #'interpolation': 'bilinear',
    }

    img = ax.imshow(data.rho.T, **options)

else:
    options = {
        'cmap': 'jet',
        'levels': 20,
        'extent' : [data.xmin, data.xmax,
                    data.ymin, data.ymax],
    }

    img = ax.contour(data.rho.T, **options)


plt.colorbar(img, cax=ax_cb)
ax_cb.yaxis.tick_right()

fig.tight_layout()

plt.show()
