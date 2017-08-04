import numpy as np
import matplotlib.pyplot as plt

# Non-Contacts
#eps = 1.0
# C C
#sig = 0.34
# N N
#sig = 0.325
# N C
#sig = (0.325+0.34)/2.
# O O
#sig = 0.296
# N O
#sig = (0.325+0.296)/2.
# O C
#sig = (0.34+0.296)/2.
# zero
#sig = 0.

# Contacts
#eps = 3.8
#sig = 0.5

#sig6 = sig**6
#sig12 = sig**12

rmin = 0.1
rcut = 0.1
rmax = 2.0
dr = 0.001
Nbins = int( (rmax-rmin)/dr ) + 2
print Nbins
r = np.linspace(rmin,rmax,Nbins)

# WCA force
#pot = [ 4.*eps*( -sig6*(x**-6) + sig12*(x**-12) ) if x < sig else 0 for x in r ]
#force = 4.*eps*( -6*sig6*(r**-7) + 12*sig12*(r**-13) )
# LJ force
#pot = eps*( -2.*sig6*(r**-6) + sig12*(r**-12) ) 
#force = eps*( -2.*6*sig6*(r**-7) + 12*sig12*(r**-13) )
# just 1/r6
pot_r6 = np.array( [ -(x**-6) if x > rcut else -(rcut**-6) for x in r ] )
force_r6 = np.array( [ -6*(x**-7) if x > rcut else -6*(rcut**-7) for x in r ] )
# just 1/r12
pot_r12 = np.array( [ (x**-12) if x > rcut else (rcut**-12) for x in r ] )
force_r12 = np.array( [ 12*(x**-13) if x > rcut else 12*(rcut**-13) for x in r ] )
#
pot_C = [ 0. for x in r ]
force_C = [ 0. for x in r ]

pot_r6 = -1.*pot_r6
force_r6 = -1.*force_r6
#
#sig = 0.3 # apply a default sig to make the force values reasonable
#sig6 = sig**6
#sig12 = sig**12
pot_r6 = (10**-3)*pot_r6
force_r6 = (10**-3)*force_r6
pot_r12 = (10**-6)*pot_r12
force_r12 = (10**-6)*force_r12

np.savetxt('pot_r6.dat', np.vstack([r, pot_r6]).T, fmt='%0.8f')
np.savetxt('force_r6.dat', np.vstack([r, force_r6]).T, fmt='%0.8f')
np.savetxt('pot_r12.dat', np.vstack([r, pot_r12]).T, fmt='%0.8f')
np.savetxt('force_r12.dat', np.vstack([r, force_r12]).T, fmt='%0.8f')
#np.savetxt('tablep.xvg', np.vstack([r, pot_C, force_C, pot_r6, force_r6, pot_r12, force_r12]).T, fmt='%0.8f')



