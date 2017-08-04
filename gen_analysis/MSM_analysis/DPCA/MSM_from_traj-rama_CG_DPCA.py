
# coding: utf-8

# CG MSM from "primitive clustering" along rama-plot
# ====

# In[1]:

import pyemma
pyemma.__version__


# In[2]:

import os
import numpy as np
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


# Read in the dtrajs
# ------

# In[4]:

Nprune = 1 # only use a portion of the data for clustering
np.save('Nprune',Nprune)

indir_1='/data/isilon/rudzinski/cluster_tmp/ALA7/ALA7_hybrid_AMBER_Go/wDB-HP_inter/Find_Tfold/NC_CA/2016_10_11/epsNC-10/epsNC-10_epsdb-0.3epsNC_epshp-0.4epsNC/T-255/'
indir_2=indir_1[:-1]+'_less-samp/'

dtraj_rama_1 = np.load(indir_1+'MSM_analysis/BMSM/dtraj/traj_rama.npy')
dtraj_rama_2 = np.load(indir_2+'MSM_analysis/BMSM/dtraj/traj_rama.npy')

#(5, 20, 400001, 2)
#(5, 20, 200001, 2)

Aconv = np.pi/180.
dtraj_dih = []
for traj in range(dtraj_rama_1.shape[1]):
    dtraj_dih.append([])
    for dih in range(dtraj_rama_1.shape[3]):
        for res in range(dtraj_rama_1.shape[0]):
            if ( (res==0) and (dih==0) ):
                dtraj_dih[traj].append(np.cos(Aconv*dtraj_rama_1[res,traj,:,dih]))
                dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.sin(Aconv*dtraj_rama_1[res,traj,:,dih])) )
            else:
                dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.cos(Aconv*dtraj_rama_1[res,traj,:,dih])) )
                dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.sin(Aconv*dtraj_rama_1[res,traj,:,dih])) )

for traj in range(dtraj_rama_2.shape[1]):
    traj_shift = traj + dtraj_rama_1.shape[1]
    dtraj_dih.append([])
    for dih in range(dtraj_rama_2.shape[3]):
        for res in range(dtraj_rama_2.shape[0]):
            if ( (res==0) and (dih==0) ):
                dtraj_dih[traj_shift].append(np.cos(Aconv*dtraj_rama_2[res,traj,:,dih]))
                dtraj_dih[traj_shift] = np.vstack( (dtraj_dih[traj_shift],np.sin(Aconv*dtraj_rama_2[res,traj,:,dih])) )
            else:
                dtraj_dih[traj_shift] = np.vstack( (dtraj_dih[traj_shift],np.cos(Aconv*dtraj_rama_2[res,traj,:,dih])) )
                dtraj_dih[traj_shift] = np.vstack( (dtraj_dih[traj_shift],np.sin(Aconv*dtraj_rama_2[res,traj,:,dih])) )

# In[6]:

for traj in range(dtraj_rama_1.shape[1]):
    dtraj_dih[traj] = dtraj_dih[traj].T
for traj in range(dtraj_rama_2.shape[1]):
    traj_shift = traj + dtraj_rama_1.shape[1]
    dtraj_dih[traj_shift] = dtraj_dih[traj_shift].T


# In[7]:

dtraj_dih[0].shape


# In[8]:

pca_obj = coor.pca(dtraj_dih, dim=-1, var_cutoff=0.95,stride=1, mean=None)


# In[9]:

Y = pca_obj.get_output()


# In[10]:

# nb - this will introduce errors into the clustering but that only matter for the mpp part, let's ignore for now
from copy import deepcopy
dtraj_conc = deepcopy(Y[0][::Nprune,0:5])
for traj in range(1,dtraj_rama_1.shape[1]):
    dtraj_conc = np.vstack((dtraj_conc,Y[traj][::Nprune,0:5]))
for traj in range(1,dtraj_rama_2.shape[1]):
    traj_shift = traj + dtraj_rama_1.shape[1]
    dtraj_conc = np.vstack((dtraj_conc,Y[traj_shift][::Nprune,0:5]))

# In[14]:

dtraj_conc.shape


# In[15]:

# save the first 5 pca dimensions
np.savetxt('dtraj_pca5D.dat',dtraj_conc)


