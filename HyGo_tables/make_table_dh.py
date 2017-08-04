import numpy as np
#import matplotlib.pyplot as plt

# Non-Contacts
er = 80.
er_ref = 80.
T_ref = 300.
Cs = 0.01
Bk = 1.
kappa = 0.32*np.sqrt(Cs)
KC = 10. # 1389
qiqj = 1.

rmin = 0.1
rcut = 0.1
rmax = 4.1
dr = 0.001
Nbins = int( (rmax-rmin)/dr ) + 2
print Nbins
r = np.linspace(rmin,rmax,Nbins)

# El force
pot = KC*Bk*qiqj*np.exp(-kappa*r) / (er*r) 
force = KC*Bk*qiqj*np.exp(-kappa*r)*(1.0/r)*( (1.0/r) + kappa )  

# apply a rescaling to counteract the default parameters
force /= 0.970906

np.savetxt('pot_dh_er-'+str(er)+'.dat', np.vstack([r, pot]).T, fmt='%0.8f')
np.savetxt('force_dh_er-'+str(er)+'.dat', np.vstack([r, force]).T, fmt='%0.8f')
#np.savetxt('table_'+type1+'_'+type2+'.xvg', np.vstack([r, pot_C, force_C, pot_r6, force_r6, pot_r12, force_r12]).T, fmt='%0.8f')



