import numpy as np
import matplotlib.pyplot as pl
import numpy.random as ra

# This script just test the random SN event generation in the SUPERNOVA module
# to make sure i am indeed getting the correct number of total supernovae and
# the correct mean number of supernovae per time-step.

pl.ion()

# Number of steps over which to explode (dur / dt in the code)
nsteps_left = 11000

# Total number of supernovae to explode
num_left = 42132

sarr = []
narr = []

# Loop over (time)steps
for i in range(int(nsteps_left) - 1):

    # probability factor
    f = num_left / nsteps_left + 1

    num = 0;
    while((f * ra.random() > 1) and (num < 1000)):

        num += 1
        num_left -= 1

    nsteps_left -= 1

    sarr.append(i)
    narr.append(num)


ntot = np.sum(narr)
print('Total number of SN: ' + format(np.sum(narr)))

nav = 1.0 * ntot / sarr[-1]
print('Average number of SN per step: ' + format(nav))

hist, bins = np.histogram(narr, 50, (0, 50))
navh = sum(hist * bins[:-1]) / sum(hist)
print('Average number of SN per step via histogram: ' + format(navh))

pl.figure(figsize=(10,4))
pl.plot(sarr, narr, '*')

pl.figure(figsize=(7,5))
pl.hist(narr, 50, (0, 50))

pl.show()

