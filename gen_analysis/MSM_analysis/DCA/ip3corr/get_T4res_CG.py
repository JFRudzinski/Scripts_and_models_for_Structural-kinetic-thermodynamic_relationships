
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
cc4res = np.load('../cc4res.npy')
N_res = len(cc_full[0])
N_quads = N_res-3

# map each full state to a sequence of 4res states
state_full_to_4res = []
for state in cc_full:
    state_full_to_4res.append( np.array([ np.where( np.all( cc4res==state[quad:quad+4], axis=1 ) )[0][0] for quad in range(N_quads) ]) )
state_full_to_4res = np.array(state_full_to_4res)

# get the traj along quads of residues
dtraj_quads = [ [] for quad in range(N_quads) ]
for traj in range(len(dtraj_DPCA)):
    for quad in range(N_quads):
        dtraj_quads[quad].append( state_full_to_4res[dtraj_DPCA[traj]][:,quad] )

# calc the count matrices
Cmat_quads = []
for quad in range(N_quads):
    Cmat_quads.append( pyemma.msm.estimation.count_matrix(dtraj_quads[quad], tau, sliding=True, sparse_return=False) )

# combine the individual count matrices (i.e., no seq dep)
Cmat = deepcopy(Cmat_quads[0])
for quad in range(1,N_quads):
    Cmat += Cmat_quads[quad]

# get the mle
lcc = pyemma.msm.estimation.largest_connected_set(Cmat, directed=True)
Cmat = pyemma.msm.estimation.largest_connected_submatrix(Cmat, directed=True, lcc=lcc)
Tmle = pyemma.msm.estimation.transition_matrix(Cmat, reversible=True)

# Now save all the model
np.save('T4res',Tmle)

