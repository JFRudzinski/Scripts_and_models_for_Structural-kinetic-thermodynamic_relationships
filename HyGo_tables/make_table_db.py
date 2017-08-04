import numpy as np
import matplotlib.pyplot as plt

# params
k = 6
m = 3
n = 2
# interaction specific!
r_cm = 0.5 # CA 1-4
eps = 10
eps_ssm = 0.005*eps 
eps_db = 0.5*eps
#
#eps_ssm = 1e-8*eps
#eps_db = 1e-6*eps
#
r_ssm = r_cm + 0.3
r_db = (r_cm + r_ssm) / 2.
#
C = 4*n*(eps+eps_db) / ((r_db-r_cm)**(4*n))
B = m*eps_ssm*((r_ssm-r_db)**(2*(m-1)))
h1 = (1.-1./m)*((r_ssm-r_db)**2) / ((eps_ssm/eps_db) + 1.)
h2 = (m-1.)*((r_ssm-r_db)**(2*m)) / (1. + (eps_db/eps_ssm))
#
rmin = 0.1
rcut = 0.1
rmax = 2.0
dr = 0.001
Nbins = int( (rmax-rmin)/dr ) + 2
print Nbins
r = np.linspace(rmin,rmax,Nbins)
#
Z_r = lambda r: (r_cm/x)**k 
#Z_r = np.array( [ (r_cm/x)**k for x in r ] )
#Z_rcut = (r_cm/rcut)**k
Zp_r = lambda r: -k*(r_cm**k)/(x**(k+1))
#Zp_r = np.array( [ k*(r_cm**k)/(x**(k+1)) for x in r ] )
#Zp_rcut = k*(r_cm**k)/(rcut**(k+1))
Y_r = lambda r: (r - r_db)**2
#Y_r = np.array( [ (r - r_db)**2 for x in r ] )
#Y_rcut = (rcut - r_db)**2
Yp_r = lambda r: 2*(r - r_db)
#Yp_r = np.array( [ 2*(r - r_db) for x in r ] )
#Yp_rcut = 2*(rcut - r_db)
#
#r = r.tolist()
#
pot_r12 = np.array( [ eps*Z_r(rcut)*(Z_r(rcut)-2.) if (x<rcut) else eps*Z_r(x)*(Z_r(x)-2.) if (x<r_cm) else C*(Y_r(x)**n)*(0.5*(Y_r(x)**n)-((r_db-r_cm)**(2*n)))/(2*n) + eps_db if (x<r_db) else -B*(Y_r(x)-h1)/((Y_r(x)**m)+h2) for x in r ] )
#pot_r12 = np.array( [ eps*Z_rcut*(Z_rcut-2.) if (x<rcut) else eps*Z_r*(Z_r-2.) if (x<r_cm) else C*(Y_r**n)*((Y_r**n)-((r_db-r_cm)**(2*n)))/(2*n) + eps_db if (x<r_db) else -B*(Y_r-h1)((Y_r**m)+h2) for x in r ] )
force_r12 = np.array( [ -eps*(Z_r(rcut)*Zp_r(rcut) + (Z_r(rcut)-2.)*Zp_r(rcut)) if (x<rcut) else -eps*(Z_r(x)*Zp_r(x) + (Z_r(x)-2.)*Zp_r(x)) if (x<r_cm) else -C*( (Y_r(x)**n)*((n/2.)*(Y_r(x)**(n-1))*Yp_r(x))+n*(Y_r(x)**(n-1))*Yp_r(x)*(0.5*(Y_r(x)**n)-((r_db-r_cm)**(2*n))) )/(2*n) if (x<r_db) else B*( ((Y_r(x)**m)+h2)*Yp_r(x) - (Y_r(x)-h1)*(m*(Y_r(x)**(m-1))*Yp_r(x)) ) / (((Y_r(x)**m)+h2)**2) for x in r ] )
#force_r12 = np.array( [ -1.0*eps*(Z_r*Zp_r + (Z_r-2.)*Zp_r) if (x<r_cm) else -C*(Y_r**n)*((Y_r**n)*Yp_r/4.)+n*(Y_r**(n-1))*Yp_r*((Y_r**n)-((r_db-r_cm)**(2*n)))/(2*n) if (x<r_db) else B*( ((Y_r**m)+h2)*Yp_r - (Y_r-h1)*(m*(Y_r**(m-1))*Yp_r) ) / (((Y_r**m)+h2)**2) for x in r ] )

# apply a rescaling to counteract the default parameters
force_r12 /= 1.67751

# just 1/r6
#pot_r6 = [ 0. for x in r ]
#force_r6 = [ 0. for x in r ]
#
#sig = 0.3 # apply a default sig to make the force values reasonable
#sig6 = sig**6
#sig12 = sig**12
#pot_r6 = (10**-3)*pot_r6
#force_r6 = (10**-3)*force_r6
#pot_r12 = (10**-6)*pot_r12
#force_r12 = (10**-6)*force_r12

#np.savetxt('pot_r6.dat', np.vstack([r, pot_r6]).T, fmt='%0.8f')
#np.savetxt('force_r6.dat', np.vstack([r, force_r6]).T, fmt='%0.8f')
np.savetxt('pot_db_rcm-'+str(r_cm)+'_eps-'+str(eps)+'_epsssm-'+str(eps_ssm)+'_epsdb-'+str(eps_db)+'.dat', np.vstack([r, pot_r12]).T, fmt='%0.8f')
np.savetxt('force_db_rcm-'+str(r_cm)+'_eps-'+str(eps)+'_epsssm-'+str(eps_ssm)+'_epsdb-'+str(eps_db)+'.dat', np.vstack([r, force_r12]).T, fmt='%0.8f')
#np.savetxt('tablep.xvg', np.vstack([r, pot_C, force_C, pot_r6, force_r6, pot_r12, force_r12]).T, fmt='%0.8f')



