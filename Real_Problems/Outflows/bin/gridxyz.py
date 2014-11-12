import numpy as np

n = 64
b = 1.024
hc = 0.5*b/n

gx = np.linspace(hc, b-hc, n).astype('f8')
gy = np.linspace(-b/2.+hc, b/2.-hc, n).astype('f8')
gz = gy

gx.tofile('gridx.out')
gy.tofile('gridy.out')
gz.tofile('gridz.out')
