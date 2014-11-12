# This script just converts Ralph's log cooling table neqstd2009.T4 
# to a linear cooling table used in PLUTO

import numpy as np
import phystools as pt

# ________________________________________________________________________________________________________
# Set this 

plot = True
write = True
br_rel_correct = True
extend = True
n_extend = 100
logT_extend = 15

# ________________________________________________________________________________________________________


# The radiative cooling tale is in log space
rct = pt.RCTable('neqstd2009.T4')

# Go to linear space.
rct.data['T'] = 10**rct.data['T']
rct.data['L'] = 10**rct.data['L']

# Find normalization constant for Thermal bremsstrahlung.
rcb = pt.RCBrems()
brh = rct.data['L'][-1]/rcb.cool(rct.data['T'][-1])
rcb = pt.RCBrems(iT0=rcb.iT0*brh)

# Begin concatenating temperature array and cooling function. 
te_new = rct.data['T']
lm_new = rct.data['L']

# Correct for relativistic Bremsstrahlung if required.
if br_rel_correct:
    rcrb = pt.RCRelBrems(iT0=rcb.iT0)
    lm_new = rcrb.cool_correct(te_new, lm_new)


# Extend cooling table to 10^15 K for Bremsstrahlung
# creating a new temperature array if required
if extend:
    te_extra = np.logspace(np.log10(rct.data['T'][-1]), logT_extend, n_extend)[1:]
    te_new = np.append(rct.data['T'], te_extra)
    if br_rel_correct:
        lm_new = np.append(lm_new, rcrb.cool(te_extra))
    else:
        lm_new = np.append(lm_new, rcb.cool(te_extra))


# Write new cooling function
if write:
    f = open('cooltable.dat','w')
    np.savetxt(f, np.vstack((te_new, lm_new)).T, fmt='%-14.6e', delimiter=' ')
    f.close()

# Plot results
if plot:
    import seaborn
    import matplotlib.pyplot as pl
    pl.plot(np.log10(rct.data['T']), np.log10(rct.data['L']), label='original')
    pl.plot(np.log10(te_new), np.log10(lm_new), label='modified')
    pl.plot(np.log10(te_new), np.log10(rcb.cool(te_new)), label='th brems')
    if br_rel_correct:
        pl.plot(np.log10(te_new), np.log10(rcrb.cool(te_new)), label='rel brems')
    pl.xlabel(r'$\log\,T\,(\mathrm{K})$', size=13)
    pl.ylabel(r'$\log\,\Lambda\,(\mathrm{erg}\,\mathrm{cm}^3\,\mathrm{s}^{-1})$', size=13)
    pl.xlim(4,10)
    pl.ylim(-24,-22)
    pl.legend(loc=4)
    pl.savefig('neqstd2009.png', bbox_inches='tight')


