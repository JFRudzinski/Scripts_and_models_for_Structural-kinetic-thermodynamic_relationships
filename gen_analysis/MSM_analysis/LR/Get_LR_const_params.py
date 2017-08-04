
# coding: utf-8

# Count the number of instances of each LR weight per residue
# ====

import pyemma
pyemma.__version__


import os
#get_ipython().magic(u'pylab inline')
#matplotlib.rcParams.update({'font.size': 12})
import numpy as np

# In[4]:

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
import msmbuilder
from msmbuilder.msm.ratematrix import ContinuousTimeMSM
import scipy
from msmtools.analysis.dense.decomposition import eigenvectors, eigenvalues
import operator


# Read in the dtrajs
# ------

# In[7]:

indir = '/data/isilon/rudzinski/cluster_tmp/ALA7/ALA7_hybrid_AMBER_Go/wDB-HP_inter/Find_Tfold/NC_CA/2016_10_11/'
state_weights = np.load(indir+'state_LR_weights.npy')
state_weights_cap = np.load(indir+'state_LR_weights_cap.npy')
state_weights_adj = np.load(indir+'state_LR_weights_adj.npy')
dtrajs = np.load('../MLE/dtraj_DPCA_CG_mapped_to_rama.npy')

# 


Nres = state_weights.shape[1]
Nw = np.zeros(Nres)
Nv = np.zeros(Nres)
Ncoil = np.zeros(Nres)
Nk = len(dtrajs)*dtrajs[0].shape[0] # total number of samples
for traj in range(len(dtrajs)):
    counts = np.bincount(dtrajs[traj])
    for state in range(counts.shape[0]):
        count = counts[state] # how many frames sample this state
        # v
        n_s = np.where(state_weights[state]==1)[0] # which residues have weight v
        if ( len(n_s) != 0 ):
            Nv[n_s] += count
        # w
        n_s = np.where(state_weights[state]==2)[0] # which residues have weight w
        if ( len(n_s) != 0 ):
            Nw[n_s] += count
        # coil
        n_s = np.where(state_weights[state]==0)[0] # which residues have weight w
        if ( len(n_s) != 0 ):
            Ncoil[n_s] += count
            
# everything is backwards (i.e., N->C-terminal), flip the arrays
#Nw = Nw[::-1]
#Nv = Nv[::-1]
#Ncoil = Ncoil[::-1]

np.savez('LR_const_params',Nw=Nw,Nv=Nv,Nk=Nk)


#


Nres_cap = state_weights_cap.shape[1]
Nw_cap = np.zeros(Nres_cap)
Nv_cap = np.zeros(Nres_cap)
Nc_cap = np.zeros(Nres_cap)
Nn_cap = np.zeros(Nres_cap)
for traj in range(len(dtrajs)):
    counts = np.bincount(dtrajs[traj])
    for state in range(counts.shape[0]):
        count = counts[state] # how many frames sample this state
        # v
        n_s = np.where(state_weights_cap[state]==1)[0] # which residues have weight v
        if ( len(n_s) != 0 ):
            Nv_cap[n_s] += count
        # w
        n_s = np.where(state_weights_cap[state]==2)[0] # which residues have weight w
        if ( len(n_s) != 0 ):
            Nw_cap[n_s] += count
        # c
        n_s = np.where(state_weights_cap[state]==4)[0] # which residues have weight c
        if ( len(n_s) != 0 ):
            Nc_cap[n_s] += count
        # n
        n_s = np.where(state_weights_cap[state]==3)[0] # which residues have weight n
        if ( len(n_s) != 0 ):
            Nn_cap[n_s] += count
    
# everything is backwards (i.e., N->C-terminal), flip the arrays
#Nw_cap = Nw_cap[::-1]
#Nv_cap = Nv_cap[::-1]
#tmp_Nn_cap = Nc_cap[::-1]
#tmp_Nc_cap = Nn_cap[::-1]
#Nn_cap = tmp_Nn_cap
#Nc_cap = tmp_Nc_cap

np.savez('LR_const_params_cap',Nw=Nw_cap,Nv=Nv_cap,Nn=Nn_cap,Nc=Nc_cap,Nk=Nk)


#


Nres = state_weights.shape[1]
Nw_adj = np.zeros(Nres)
Nv_adj = np.zeros(Nres)
Nk = len(dtrajs)*dtrajs[0].shape[0] # total number of samples
for traj in range(len(dtrajs)):
    counts = np.bincount(dtrajs[traj])
    for state in range(counts.shape[0]):
        count = counts[state] # how many frames sample this state
        # v
        n_s = np.where(state_weights_adj[state]==1)[0] # which residues have weight v
        if ( len(n_s) != 0 ):
            Nv_adj[n_s] += count
        # w
        n_s = np.where(state_weights_adj[state]==2)[0] # which residues have weight w
        if ( len(n_s) != 0 ):
            Nw_adj[n_s] += count
            
# everything is backwards (i.e., N->C-terminal), flip the arrays
#Nw_adj = Nw_adj[::-1]
#Nv_adj = Nv_adj[::-1]

np.savez('LR_const_params_adj',Nw=Nw_adj,Nv=Nv_adj,Nk=Nk)


