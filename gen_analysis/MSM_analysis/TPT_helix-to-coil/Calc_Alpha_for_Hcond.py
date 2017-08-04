
# coding: utf-8

# In[ ]:

import pyemma
pyemma.__version__

import numpy as np

# In[ ]:

import os
#get_ipython().magic(u'pylab inline')
#matplotlib.rcParams.update({'font.size': 12})
#from matplotlib import rc
#rc('font', **{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)


# In[ ]:

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
import msmbuilder
from msmbuilder.msm.ratematrix import ContinuousTimeMSM
import scipy
from msmtools.analysis.dense.decomposition import eigenvectors, eigenvalues
import operator


# Get the model
# ------

TNres_CG = np.load('../MLE/T_mle.npy')

from copy import deepcopy
s = 0
d = 31
# Get the probability that a traj passes from i to j through d
Alpha_CG = np.zeros(shape=(TNres_CG.shape[0],TNres_CG.shape[1]))
for row in range(TNres_CG.shape[0]):
    for col in np.delete(np.arange(TNres_CG.shape[0]),[row]):
        # Get Tbar
        Tbar = deepcopy(TNres_CG)
        # make the j and d states sinks
        Tbar[col] *= 0.
        Tbar[col,col] = 1.
        Tbar[d] *= 0.
        Tbar[d,d] = 1.
        #
        evals, evecs = scipy.linalg.eig(Tbar)
        evals = np.real(evals)
        evecs = np.real(evecs)
        sinks_ind = np.where(np.abs(1.-evals)<1e-6)[0]
        if ( len(sinks_ind) > 2 ):
            raise ValueError('More than 2 sinks found')
        evecs = evecs[:,sinks_ind]
        if ( col != d ): # find the evec corresponding to prob of reaching d before the other sink
            d_ind = np.where(evecs[col]<1e-6)[0]
            if ( len(d_ind) != 1 ):
                raise ValueError('Could not identify the evec corresponding to the coil state')
            Alpha_CG[row,col] = evecs[row,d_ind] / np.max(evecs[:,d_ind])
        else: # there is only one sink, this value should never be used anyway
            Alpha_CG[row,col] = evecs[row] / np.max(evecs)

np.save('Alpha_unf',Alpha_CG)

s = 31
d = 0
# Get the probability that a traj passes from i to j through d
Alpha_CG = np.zeros(shape=(TNres_CG.shape[0],TNres_CG.shape[1]))
for row in range(TNres_CG.shape[0]):
    for col in np.delete(np.arange(TNres_CG.shape[0]),[row]):
        # Get Tbar
        Tbar = deepcopy(TNres_CG)
        # make the j and d states sinks
        Tbar[col] *= 0.
        Tbar[col,col] = 1.
        Tbar[d] *= 0.
        Tbar[d,d] = 1.
        #
        evals, evecs = scipy.linalg.eig(Tbar)
        evals = np.real(evals)
        evecs = np.real(evecs)
        sinks_ind = np.where(np.abs(1.-evals)<1e-6)[0]
        if ( len(sinks_ind) > 2 ):
            raise ValueError('More than 2 sinks found')
        evecs = evecs[:,sinks_ind]
        if ( col != d ): # find the evec corresponding to prob of reaching d before the other sink
            d_ind = np.where(evecs[col]<1e-6)[0]
            if ( len(d_ind) != 1 ):
                raise ValueError('Could not identify the evec corresponding to the coil state')
            Alpha_CG[row,col] = evecs[row,d_ind] / np.max(evecs[:,d_ind])
        else: # there is only one sink, this value should never be used anyway
            Alpha_CG[row,col] = evecs[row] / np.max(evecs)

np.save('Alpha_fol',Alpha_CG)


