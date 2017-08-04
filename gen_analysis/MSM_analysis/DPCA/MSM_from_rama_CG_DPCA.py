
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

indir='/data/isilon/rudzinski/cluster_tmp/ALA7/ALA7_hybrid_AMBER_Go/wDB-HP_inter/Find_Tfold/NC_CA/2016_10_11/epsNC-11/epsNC-11_epsdb-0.05epsNC_epshp-0.6epsNC/T-250/'

traj_dir_base = 'run_from_Qhel-'
Qhel_val = ['0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9']

dtraj_rama_2 = []
dtraj_rama_3 = []
dtraj_rama_4 = []
dtraj_rama_5 = []
dtraj_rama_6 = []
for i in range (len(Qhel_val)):
    for j in range(2):
        traj_dir = indir+traj_dir_base+Qhel_val[i]+'/run'+str(j+1)+'/'
        dtraj_rama_2.append( np.genfromtxt(traj_dir+'rama_ALA2.xvg')[:,0:2] )
        dtraj_rama_3.append( np.genfromtxt(traj_dir+'rama_ALA3.xvg')[:,0:2] )
        dtraj_rama_4.append( np.genfromtxt(traj_dir+'rama_ALA4.xvg')[:,0:2] )
        dtraj_rama_5.append( np.genfromtxt(traj_dir+'rama_ALA5.xvg')[:,0:2] )
        dtraj_rama_6.append( np.genfromtxt(traj_dir+'rama_ALA6.xvg')[:,0:2] )


# In[5]:

Aconv = np.pi/180.
dtraj_dih = []
for traj in range(len(dtraj_rama_2)):
    dtraj_dih.append([])
    dtraj_dih[traj].append(np.cos(Aconv*dtraj_rama_2[traj][:,0]))
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.sin(Aconv*dtraj_rama_2[traj][:,0])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.cos(Aconv*dtraj_rama_3[traj][:,0])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.sin(Aconv*dtraj_rama_3[traj][:,0])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.cos(Aconv*dtraj_rama_4[traj][:,0])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.sin(Aconv*dtraj_rama_4[traj][:,0])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.cos(Aconv*dtraj_rama_5[traj][:,0])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.sin(Aconv*dtraj_rama_5[traj][:,0])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.cos(Aconv*dtraj_rama_6[traj][:,0])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.sin(Aconv*dtraj_rama_6[traj][:,0])) )
    #
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.cos(Aconv*dtraj_rama_2[traj][:,1])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.sin(Aconv*dtraj_rama_2[traj][:,1])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.cos(Aconv*dtraj_rama_3[traj][:,1])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.sin(Aconv*dtraj_rama_3[traj][:,1])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.cos(Aconv*dtraj_rama_4[traj][:,1])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.sin(Aconv*dtraj_rama_4[traj][:,1])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.cos(Aconv*dtraj_rama_5[traj][:,1])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.sin(Aconv*dtraj_rama_5[traj][:,1])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.cos(Aconv*dtraj_rama_6[traj][:,1])) )
    dtraj_dih[traj] = np.vstack( (dtraj_dih[traj],np.sin(Aconv*dtraj_rama_6[traj][:,1])) )


# In[6]:

for traj in range(len(dtraj_rama_2)):
    dtraj_dih[traj] = dtraj_dih[traj].T


# In[7]:

dtraj_dih[0].shape


# In[8]:

pca_obj = coor.pca(dtraj_dih, dim=-1, var_cutoff=0.95,stride=1, mean=None)


# In[9]:

Y = pca_obj.get_output()


# In[10]:

# nb - this will introduce errors into the clustering but that only matter for the mpp part, let's ignore for now
from copy import deepcopy
dtraj_conc = deepcopy(Y[0][:,0:5])
for traj in range(1,len(dtraj_rama_2)):
    dtraj_conc = np.vstack((dtraj_conc,Y[traj][:,0:5]))


# In[14]:

dtraj_conc.shape


# In[15]:

# save the first 5 pca dimensions
np.savetxt('dtraj_pca5D.dat',dtraj_conc)


