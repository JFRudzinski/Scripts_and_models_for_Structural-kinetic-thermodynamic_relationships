import numpy as np
import matplotlib.pyplot as plt

# Non-Contacts
eps = 1.0
# C C
#type1 = 'C' 
#type2 = 'C'
#sig = 0.34
#sig /= (2.**(1/6.))
# N N
#type1 = 'N' 
#type2 = 'N'
#sig = 0.325
#sig /= (2.**(1/6.))
# N C
#type1 = 'N'
#type2 = 'C'
#sig = (0.325+0.34)/2.
#sig /= (2.**(1/6.))
# O O
#type1 = 'O'
#type2 = 'O'
#sig = 0.296
#sig /= (2.**(1/6.))
# N O
#type1 = 'N'
#type2 = 'O'
#sig = (0.325+0.296)/2.
#sig /= (2.**(1/6.))
# O C
type1 = 'C'
type2 = 'O'
sig = (0.34+0.296)/2.
sig /= (2.**(1/6.))
# zero
#sig = 0.

# Contacts
#eps = 3.8
#sig = 0.5

sig6 = sig**6
sig12 = sig**12

rmin = 0.1
rcut = 0.1
rmax = 2.0
dr = 0.001
Nbins = int( (rmax-rmin)/dr ) + 2
print Nbins
r = np.linspace(rmin,rmax,Nbins)


# WCA force
pot_r12_cut = 4.*eps*( -sig6*(rcut**-6) + sig12*(rcut**-12) + 1. )
pot_r12 = [ 4.*eps*( -sig6*(x**-6) + sig12*(x**-12) + 1/4.) if (x < (2.**(1/6.)*sig) and x >= rcut) else pot_r12_cut if (x < rcut) else 0 for x in r ]
force_r12_cut = 4.*eps*( -6*sig6*(rcut**-7) + 12*sig12*(rcut**-13) )
force_r12 = [ 4.*eps*( -6*sig6*(x**-7) + 12*sig12*(x**-13) ) if (x < (2.**(1/6.)*sig) and x >= rcut) else force_r12_cut if (x < rcut)  else 0 for x in r ]
pot_r6 = [ 0. for x in r ]
force_r6 = [ 0. for x in r ]
# LJ force
#pot = eps*( -2.*sig6*(r**-6) + sig12*(r**-12) ) 
#force = eps*( -2.*6*sig6*(r**-7) + 12*sig12*(r**-13) )
# just 1/r6
#pot_r6 = [ -(x**-6) if x > rcut else -(rcut**-6) for x in r ]
#force_r6 = [ -6*(x**-7) if x > rcut else -6*(rcut**-7) for x in r ]
# just 1/r12
#pot_r12 = [ (x**-12) if x > rcut else (rcut**-12) for x in r ]
#force_r12 = [ 12*(x**-13) if x > rcut else 12*(rcut**-13) for x in r ]
#
pot_C = [ 0. for x in r ]
force_C = [ 0. for x in r ]

np.savetxt('pot_'+type1+'_'+type2+'.dat', np.vstack([r, pot_r12]).T, fmt='%0.8f')
np.savetxt('force_'+type1+'_'+type2+'.dat', np.vstack([r, force_r12]).T, fmt='%0.8f')
#np.savetxt('table_'+type1+'_'+type2+'.xvg', np.vstack([r, pot_C, force_C, pot_r6, force_r6, pot_r12, force_r12]).T, fmt='%0.8f')



