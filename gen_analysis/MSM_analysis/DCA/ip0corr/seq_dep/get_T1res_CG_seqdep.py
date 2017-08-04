
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
N_res = len(cc_full[0])

# get the traj along individual residues
dtraj_res = [ [] for res in range(N_res) ]
for traj in range(len(dtraj_DPCA)):
    for res in range(N_res):
        dtraj_res[res].append( cc_full[dtraj_DPCA[traj]][:,res] )

# calc the count matrices
Cmat_res = []
for res in range(N_res):
    Cmat_res.append( pyemma.msm.estimation.count_matrix(dtraj_res[res], tau, sliding=True, sparse_return=False) )

# get the mle
Tmle = []
for res in range(N_res):
    lcc = pyemma.msm.estimation.largest_connected_set(Cmat_res[res], directed=True)
    Cmat_res[res] = pyemma.msm.estimation.largest_connected_submatrix(Cmat_res[res], directed=True, lcc=lcc)
    Tmle.append( pyemma.msm.estimation.transition_matrix(Cmat_res[res], reversible=True) )

# Now save all the model
np.save('T1res_seqdep',Tmle)

