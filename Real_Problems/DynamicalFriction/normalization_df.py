import numpy as np
import physconst as pc
import norm



# Set these
mach = 1
mbh = 1.e6
pres = 1.e6

gamma = 1.6

G = pc.newton
m = mbh * pc.msun
p = pres * pc.kboltz

cs = pow(4 * gamma * G*G*G * m*m * p * np.pi / 3, 0.125)
rdf = pc.newton * m / (cs * cs)
dens = m / (4 * np.pi / 3.) / (rdf * rdf * rdf)

print('Normalizations:')
print('cs = ' + str(cs / 1.e5) + ' km / s')
print('dens = ' + str(dens / ( pc.msun / pc.pc**3)) + ' msun / pc^3')
print('rdf = ' + str(rdf / pc.pc) + ' pc')

