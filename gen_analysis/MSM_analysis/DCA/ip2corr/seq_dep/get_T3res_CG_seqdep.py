
# coding: utf-8

# AA MSM from Qhel directly via rama-plot
# -------
# with nearest neighbor correlations
# ====

# In[1]:

import pyemma
pyemma.__version__

import numpy as np

# In[2]:

import os
#get_ipython().magic(u'pylab inline')
#matplotlib.rcParams.update({'font.size': 12})


# In[3]:

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
import msmbuilder
from msmbuilder.msm.ratematrix import ContinuousTimeMSM
import scipy
from msmtools.analysis.dense.decomposition import eigenvectors, eigenvalues
import operator
from copy import deepcopy

# Read in the dtrajs
# ------

# get the indir from the current dir
cdir = os.getcwd()
flag = False
wdir = 'MSM_analysis'
dir_len = len(wdir)
ctr = 0
while (not flag):
    if ( cdir[ctr:dir_len+ctr] == wdir ):
        flag = True
        dir_ind = ctr
    ctr += 1

indir = cdir[:ctr+dir_len]


dtraj_DPCA = np.load(indir+'MLE/dtraj_DPCA_CG_mapped_to_rama.npy')
tau = np.genfromtxt(indir+'MLE/tau_CG.dat').astype(int)
cc_full = np.load(indir+'BMSM/dtraj/cc_full.npy')
cc3res = np.load('../../cc3res.npy')
N_res = len(cc_full[0])
N_trips = N_res-2

# map each full state to a sequence of 3res states
state_full_to_3res = []
for state in cc_full:
    state_full_to_3res.append( np.array([ np.where( np.all( cc3res==state[trip:trip+3], axis=1 ) )[0][0] for trip in range(N_trips) ]) )
state_full_to_3res = np.array(state_full_to_3res)

# get the traj along trips of residues
dtraj_trips = [ [] for trip in range(N_trips) ]
for traj in range(len(dtraj_DPCA)):
    for trip in range(N_trips):
        dtraj_trips[trip].append( state_full_to_3res[dtraj_DPCA[traj]][:,trip] )

# calc the count matrices
Cmat_trips = []
for trip in range(N_trips):
    Cmat_trips.append( pyemma.msm.estimation.count_matrix(dtraj_trips[trip], tau, sliding=True, sparse_return=False) )

# get the mle
Tmle = []
for trip in range(N_trips):
    lcc = pyemma.msm.estimation.largest_connected_set(Cmat_trips[trip], directed=True)
    Cmat_trips[trip] = pyemma.msm.estimation.largest_connected_submatrix(Cmat_trips[trip], directed=True, lcc=lcc)
    Tmle.append( pyemma.msm.estimation.transition_matrix(Cmat_trips[trip], reversible=True) )

# Now save all the model
np.save('T3res_seqdep',Tmle)

