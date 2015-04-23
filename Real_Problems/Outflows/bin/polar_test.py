# This script test the CART2POL2 macro defined in macros_usr.h
# The result should be a circular pattern of increasing values of the
# angle from 0 to 2pi. The script tests whether the macro handles positive
# and negative values of x1 and x2 correctly.

import numpy as np
import matplotlib.pyplot as pl

CONST_PI = np.pi

pl.ion()

def HS(x2):
    return np.where(np.greater(x2, 0), 1, 0)

def SGN(x2):
    return np.where(np.greater(x2, 0), 1, np.where(np.less(x2, 0), -1, 0))

def CART2POL1(x1, x2, x3):
    return np.sqrt((x1)*(x1) + (x2)*(x2))

def CART2POL2(x1, x2, x3):
    return 2*CONST_PI*HS(-(x2)) + (SGN(x2))*np.arccos((x1)/(CART2POL1(x1, x2, x3)))

x1 = np.arange(-1, 1.1, 0.1)

pl.pcolormesh(x1, x1, CART2POL2(x1[:,None], x1[None,:], 0))
pl.colorbar()

