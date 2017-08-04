import numpy as np
#import matplotlib.pyplot as plt

# Non-Contacts
eps = 10
sig = 0.49592

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
pot = 4.*eps*( -sig6*(r**-6) + sig12*(r**-12) ) 
force = 4.*eps*( -6*sig6*(r**-7) + 12*sig12*(r**-13) )

# apply a rescaling to counteract the default parameters
force /= 1.67751

np.savetxt('pot_hp_eps-'+str(eps)+'.dat', np.vstack([r, pot]).T, fmt='%0.8f')
np.savetxt('force_hp_eps-'+str(eps)+'.dat', np.vstack([r, force]).T, fmt='%0.8f')
#np.savetxt('table_'+type1+'_'+type2+'.xvg', np.vstack([r, pot_C, force_C, pot_r6, force_r6, pot_r12, force_r12]).T, fmt='%0.8f')



