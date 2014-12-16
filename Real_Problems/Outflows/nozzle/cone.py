import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

pl.ion()

def conez(x, y, alpha=15, r=1.):
    return np.sqrt((x*x + y*y)/np.tan(np.deg2rad(alpha))**2)


x = np.linspace(-3,3.,61)
x2 = x[None,:]

y = np.linspace(-3,3.,61)
y2 = y[:,None]

one2 = np.ones((y.size, x.size))

x2 = one2*x2
y2 = one2*y2

fig = pl.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(x2, y2, conez(x2, y2))

