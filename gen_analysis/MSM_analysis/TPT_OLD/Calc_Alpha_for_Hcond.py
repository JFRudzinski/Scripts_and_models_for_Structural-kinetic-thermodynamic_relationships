
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


# Let's compute the conditional path entropy for each path and then perform a weighted sum according to the fractional flux
from copy import deepcopy
s = 0
d = 31
# we are going to need the probability that a traj passes from i to j through d
Alpha_CG = np.zeros(shape=(TNres_CG.shape[0],TNres_CG.shape[1]))
for row in range(TNres_CG.shape[0]):
    for col in np.delete(np.arange(TNres_CG.shape[0]),[row]):
        mle_Nres_CG = pyemma.msm.markov_model(TNres_CG)
        tpt_tmp = msm.tpt(mle_Nres_CG,[row],[col])
        try: # some of the node pairs have too low connections for this analysis??
            paths_tmp = np.array(tpt_tmp.pathways())
            # grab only paths that pass through d
            subpaths = np.where( np.array([d in path for path in paths_tmp[0]]) == True )
            Alpha_CG[row,col] = np.sum(paths_tmp[1][subpaths]) / np.sum(paths_tmp[1])
        except:
            Alpha_CG[row,col] = 0.

np.save('Alpha_unf',Alpha_CG)

s = 31
d = 0
# we are going to need the probability that a traj passes from i to j through d
Alpha_CG = np.zeros(shape=(TNres_CG.shape[0],TNres_CG.shape[1]))
for row in range(TNres_CG.shape[0]):
    for col in np.delete(np.arange(TNres_CG.shape[0]),[row]):
        mle_Nres_CG = pyemma.msm.markov_model(TNres_CG)
        tpt_tmp = msm.tpt(mle_Nres_CG,[row],[col])
        try: # some of the node pairs have too low connections for this analysis??
            paths_tmp = np.array(tpt_tmp.pathways())
            # grab only paths that pass through d
            subpaths = np.where( np.array([d in path for path in paths_tmp[0]]) == True )
            Alpha_CG[row,col] = np.sum(paths_tmp[1][subpaths]) / np.sum(paths_tmp[1])
        except:
            Alpha_CG[row,col] = 0.

np.save('Alpha_fol',Alpha_CG)


