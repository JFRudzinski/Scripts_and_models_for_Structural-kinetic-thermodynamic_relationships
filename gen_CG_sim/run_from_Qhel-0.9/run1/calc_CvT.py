import numpy as np

T = 267
kB = 0.00831451

time, Etot = np.genfromtxt('energy.xvg', skip_header=19).T

Eavg = np.sum(Etot) / Etot.shape[0]
Eavg2 = np.sum(Etot**2 / Etot.shape[0])

CvT = (Eavg2 - Eavg**2) / (kB*T**2)

np.savetxt('CvT.dat',np.c_[kB*T,CvT])

# and get the E dist
E_hist, bin_edges = np.histogram( Etot,bins=100,normed=1 )
bins = 0.5*(bin_edges[1:] + bin_edges[:-1])
np.savetxt('E_dist.dat',np.c_[bins,E_hist])
