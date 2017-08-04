
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
from copy import deepcopy

# Get the model
TNres_CG = np.load('../MLE/T_mle.npy')
# Get the probability that a traj passes from i to j through d
Alpha_CG = np.load('Alpha_fol.npy')

# Let's compute the conditional path entropy for each of the dominant reaction pathways
s = 31
d = 0

# get the reaction pathways from s to d
mle_Nres_CG = pyemma.msm.markov_model(TNres_CG)
tpt_mle_CG = msm.tpt(mle_Nres_CG,[s],[d])
paths_CG_unf = np.array(tpt_mle_CG.pathways())

# compute the conditional path entropy according to Kafsi et al. IEEE Inf Theory (2013)
H_path_cond_CG = []
for path in paths_CG_unf[0]:
    tmp = 1.
    H_tmp = 0.
    for node in range(len(path)-2):
        # Get Tp
        Tp = deepcopy(TNres_CG)
        for row in range(TNres_CG.shape[0]):
            if ( (row==path[node+1]) or (row==d) ):
                Tp[row] *= 0.
                Tp[row,row] = 1.
                continue
            elif ( abs(Alpha_CG[row,path[node+1]]-1.)<1e-6 ):
                continue
            
            for col in range(TNres_CG.shape[1]):
                Tp[row,col] *= (1.-Alpha_CG[col,path[node+1]]) / (1.-Alpha_CG[row,path[node+1]])
        # calc H(Tp)k,k+1
        Qd = np.delete(TNres_CG,[path[node+1]],axis=0) # nb - I thought this should be Tp, but then you will always have a singular matrix...confusing notation!
        Qd = np.delete(Qd,[path[node+1]],axis=1)
        I_N = np.identity(TNres_CG.shape[0]-1)
        M_F = np.linalg.inv(I_N-Qd)
        for ind,node2 in enumerate( np.delete(np.arange(TNres_CG.shape[0]),[path[node+1]]) ):
            row2 = Tp[node2]
            row2 = np.delete(row2,[node2])
            H_row = -row2*np.ma.log(row2)
            H_row = np.ma.MaskedArray(H_row,fill_value=0)
            H_row = np.ma.filled(H_row)
            m_ind = deepcopy(ind)
            if ( m_ind > path[node+1] ):
                m_ind -= 1
            s_ind = deepcopy(path[node])
            if ( s_ind > path[node+1] ):
                s_ind -= 1
            H_tmp += M_F[s_ind,m_ind]*np.sum(H_row)
    # calc H(T)ld
    l = len(path)-2
    Qd = np.delete(TNres_CG,[d],axis=0)
    Qd = np.delete(Qd,[d],axis=1)
    I_N = np.identity(TNres_CG.shape[0]-1)
    M_F = np.linalg.inv(I_N-Qd)
    for ind,node2 in enumerate( np.delete(np.arange(TNres_CG.shape[0]),[d]) ):
        row2 = TNres_CG[node2]
        row2 = np.delete(row2,[node2])
        H_row = -row2*np.ma.log(row2)
        H_row = np.ma.MaskedArray(H_row,fill_value=0)
        H_row = np.ma.filled(H_row)
        m_ind = deepcopy(ind)
        if ( m_ind > d ):
            m_ind -= 1
        s_ind = deepcopy(path[l])
        if ( s_ind > d ):
            s_ind -= 1
        H_tmp += M_F[s_ind,m_ind]*np.sum(H_row)
    H_path_cond_CG.append(H_tmp)

np.save('H_path_cond_CG_fol.npy', H_path_cond_CG)



# repeat the calculation for single intermediate states; i.e., H_{sd|u} for each u
# this is the entropy of reactive trajectories given that they pass through u
H_path_cond_state_CG = np.zeros(TNres_CG.shape[0])
for state in range(1,TNres_CG.shape[0]-1):
    path = np.array([s,state,d])
    tmp = 1.
    H_tmp = 0.
    for node in range(len(path)-2):
        # Get Tp
        Tp = deepcopy(TNres_CG)
        for row in range(TNres_CG.shape[0]):
            if ( (row==path[node+1]) or (row==d) ):
                Tp[row] *= 0.
                Tp[row,row] = 1.
                continue
            elif ( abs(Alpha_CG[row,path[node+1]]-1.)<1e-6 ):
                continue

            for col in range(TNres_CG.shape[1]):
                Tp[row,col] *= (1.-Alpha_CG[col,path[node+1]]) / (1.-Alpha_CG[row,path[node+1]])
        # calc H(Tp)k,k+1
        Qd = np.delete(TNres_CG,[path[node+1]],axis=0)
        Qd = np.delete(Qd,[path[node+1]],axis=1)
        I_N = np.identity(TNres_CG.shape[0]-1)
        M_F = np.linalg.inv(I_N-Qd)
        for ind,node2 in enumerate( np.delete(np.arange(TNres_CG.shape[0]),[path[node+1]]) ):
            row2 = Tp[node2]
            row2 = np.delete(row2,[node2])
            H_row = -row2*np.ma.log(row2)
            H_row = np.ma.MaskedArray(H_row,fill_value=0)
            H_row = np.ma.filled(H_row)
            m_ind = deepcopy(ind)
            if ( m_ind > path[node+1] ):
                m_ind -= 1
            s_ind = deepcopy(path[node])
            if ( s_ind > path[node+1] ):
                s_ind -= 1
            H_tmp += M_F[s_ind,m_ind]*np.sum(H_row)
    # calc H(T)ld
    l = len(path)-2
    Qd = np.delete(TNres_CG,[d],axis=0)
    Qd = np.delete(Qd,[d],axis=1)
    I_N = np.identity(TNres_CG.shape[0]-1)
    M_F = np.linalg.inv(I_N-Qd)
    for ind,node2 in enumerate( np.delete(np.arange(TNres_CG.shape[0]),[d]) ):
        row2 = TNres_CG[node2]
        row2 = np.delete(row2,[node2])
        H_row = -row2*np.ma.log(row2)
        H_row = np.ma.MaskedArray(H_row,fill_value=0)
        H_row = np.ma.filled(H_row)
        m_ind = deepcopy(ind)
        if ( m_ind > d ):
            m_ind -= 1
        s_ind = deepcopy(path[l])
        if ( s_ind > d ):
            s_ind -= 1
        H_tmp += M_F[s_ind,m_ind]*np.sum(H_row)
    H_path_cond_state_CG[state] = 1.*H_tmp

np.save('H_path_cond_state_CG_fol.npy', H_path_cond_state_CG)


