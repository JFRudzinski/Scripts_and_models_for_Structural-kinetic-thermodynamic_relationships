
# coding: utf-8

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


# In[4]:

# get the DPCA dtrajs, mapped to the rama states
TN = np.load('TN_opt.npy')
cc_full = np.load('../cc_full.npy')


# In[5]:

# get the stat dist of the model
mle_Nres = pyemma.msm.markov_model(TN)
mu_Nres = mle_Nres.eigenvectors_left(k=1)[0]


# In[6]:

# do the mapping
N_Qhel = 6
TQhel = np.zeros(shape=(N_Qhel,N_Qhel))
PI = np.zeros(N_Qhel)
for i in range(cc_full.shape[0]):
    Xi = np.sum(cc_full[i]).astype(int)
    PI[Xi] += mu_Nres[i]
    for j in range(cc_full.shape[0]):
        Xj = np.sum(cc_full[j]).astype(int)
        TQhel[Xi,Xj] += mu_Nres[i]*TN[i][j]
for I in range(N_Qhel):
    TQhel[I] /= PI[I]
# need to divide by the normalization factor here as usual!


# In[7]:

# choose tau and make the model
np.save('TN_opt_Qhel',TQhel)


