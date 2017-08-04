
# coding: utf-8

# Collect the data from LR sampling
# ====

# In[1]:

import pyemma
pyemma.__version__


# In[2]:

import os
#get_ipython().magic(u'pylab inline')
#matplotlib.rcParams.update({'font.size': 12})
import numpy as np

# In[3]:

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
import msmbuilder
from msmbuilder.msm.ratematrix import ContinuousTimeMSM
import scipy
from msmtools.analysis.dense.decomposition import eigenvectors, eigenvalues
import operator
import scipy


# Read in the data
# ------

# The LR sampling output data
Nproc = 16
lnw = []
lnv = []
L = []
for i in range(Nproc):
    data = np.load('output_variables_'+str(i)+'.npz')
    lnw.append(data['lnw'])
    lnv.append(data['lnv'])
    L.append(data['L'])
lnw = np.array(lnw)
lnv = np.array(lnv)
L = np.array(L)

# The LR sampling input data
LR_const_params_adj = np.load('../LR_const_params_adj.npz')
Nk = LR_const_params_adj['Nk']
Nv_adj = LR_const_params_adj['Nv']
Nw_adj = LR_const_params_adj['Nw']
#
LR_const_params = np.load('../LR_const_params.npz')
Nv = LR_const_params['Nv']
Nw = LR_const_params['Nw']

# The eqm prop weights
indir = '/data/isilon/rudzinski/cluster_tmp/ALA7/ALA7_hybrid_AMBER_Go/wDB-HP_inter/Find_Tfold/NC_CA/2016_10_11/'
state_weights_eqm_prop = np.load(indir+'state_weights_eqm_prop.npz')
Nhi_weight = state_weights_eqm_prop['Nhi_weight']
Nh_weight = state_weights_eqm_prop['Nh_weight']
Nhb_weight = state_weights_eqm_prop['Nhb_weight']
Nl_weight = state_weights_eqm_prop['Nl_weight']
l_weight = state_weights_eqm_prop['l_weight']
fh_weight = state_weights_eqm_prop['fh_weight']
Ns_weight = state_weights_eqm_prop['Ns_weight']
hi_weight = state_weights_eqm_prop['hi_weight']

cc_full = np.load('../../BMSM/dtraj/cc_full.npy')
state_LR_weights = np.load(indir+'state_LR_weights.npy')

# The MSM
Tmle = np.load('../../MLE/T_mle.npy')

proc_max, samp_max = np.where(L==np.max(L))

def wv_to_M(N_res,lnw,lnv):
    M = np.zeros(shape=(N_res,3,3))
    for res in range(N_res):
        w = np.exp(lnw[res])
        v = np.exp(lnv[res])
        M[res,0,0] = w
        M[res,0,1] = v
        M[res,1,2] = 1.
        M[res,2,0] = v
        M[res,2,1] = v
        M[res,2,2] = 1.
    return M

def M_to_Z(N_res,M):
    vec0 = np.array([0.,0.,1.])
    vecf = np.array([0.,1.,1.])
    Z = np.dot(vec0,M[0])
    for res in range(1,N_res):
        Z = np.dot(Z,M[res])
    Z = np.dot(Z,vecf)
    return Z

def calc_h_i(N_res,M,i,w):
    from copy import deepcopy
    M_i = np.zeros(shape=M[i].shape)
    M_i[0,0] = deepcopy(M[i][0,0])/w
    vec0 = np.array([0.,0.,1.])
    vecf = np.array([0.,1.,1.])
    Z = M_to_Z(N_res,M)
    if ( i == 0 ):
        dZdw = np.dot(vec0,M_i)    
    else:
        dZdw = np.dot(vec0,M[0])
    for res in range(1,N_res):
        if (i==res):
            dZdw = np.dot(dZdw,M_i)
        else:
            dZdw = np.dot(dZdw,M[res])
    dZdw = np.dot(dZdw,vecf)
    
    return (w/Z)*dZdw

def get_h(N_res,M,lnw):
    hi = []
    for res in range(N_res):
        w = np.exp(lnw[res])
        hi.append(calc_h_i(N_res,M,res,w))
    return np.array(hi)

def logL(N_res,Nw,Nv,Nk,lnw,lnv):
    L = 0.
    M = wv_to_M(N_res,lnw,lnv)
    Z = M_to_Z(N_res,M)
    for res in range(N_res):
        L += Nw[res]*lnw[res]
        L += Nv[res]*lnv[res]
        L -= Nk*np.log(Z)
    return L, Z

def calc_Ns_i(N_res,M,i,w,v):
    # we are calculating d^2 lnZ / (dlnw_i-1 dlnvi)
    from copy import deepcopy
    M_im1 = np.zeros(shape=M[i-1].shape)
    M_im1[0,0] = deepcopy(M[i-1][0,0])/w
    M_i = np.zeros(shape=M[i].shape)
    M_i[0,1] = deepcopy(M[i][0,1])/v
    M_i[2,0] = deepcopy(M[i][2,0])/v
    M_i[2,1] = deepcopy(M[i][2,1])/v
    vec0 = np.array([0.,0.,1.])
    vecf = np.array([0.,1.,1.])
    Z = M_to_Z(N_res,M)
    # dZ/dlnw_i-1
    if ( i-1 == 0 ):
        dZdlnw = np.dot(vec0,M_im1)    
    else:
        dZdlnw = np.dot(vec0,M[0])
    for res in range(1,N_res):
        if (i-1==res):
            dZdlnw = np.dot(dZdlnw,M_im1)
        else:
            dZdlnw = np.dot(dZdlnw,M[res])
    dZdlnw = np.dot(dZdlnw,vecf)
    dZdlnw = w*dZdlnw
    # dZ/dlnv_i
    if ( i == 0 ):
        dZdlnv = np.dot(vec0,M_i)    
    else:
        dZdlnv = np.dot(vec0,M[0])
    for res in range(1,N_res):
        if (i==res):
            dZdlnv = np.dot(dZdlnv,M_i)
        else:
            dZdlnv = np.dot(dZdlnv,M[res])
    dZdlnv = np.dot(dZdlnv,vecf)
    dZdlnv = v*dZdlnv
    # d2Z/dlnwlnv
    if ( i-1 == 0 ):
        d2Z = np.dot(vec0,M_im1)    
    else:
        d2Z = np.dot(vec0,M[0])
    for res in range(1,N_res):
        if (i-1==res):
            d2Z = np.dot(d2Z,M_im1)
        elif (i==res):
            d2Z = np.dot(d2Z,M_i)
        else:
            d2Z = np.dot(d2Z,M[res])
    d2Z = np.dot(d2Z,vecf)
    d2Z = w*v*d2Z
    
    return (1./(Z**2))*( (Z*d2Z) - (dZdlnw*dZdlnv) )

def get_Ns(N_res,M,lnw,lnv):
    Ns = 0.
    for res in range(2,N_res):
        w = np.exp(lnw[res-1])
        v = np.exp(lnv[res])
        Ns += calc_Ns_i(N_res,M,res,w,v)
    return Ns

def calc_Nl_i(N_res,M,i,v):
    from copy import deepcopy
    M_i = np.zeros(shape=M[i].shape)
    M_i[2,1] = deepcopy(M[i][2,1])/v
    vec0 = np.array([0.,0.,1.])
    vecf = np.array([0.,1.,1.])
    Z = M_to_Z(N_res,M)
    if ( i == 0 ):
        dZdv = np.dot(vec0,M_i)    
    else:
        dZdv = np.dot(vec0,M[0])
    for res in range(1,N_res):
        if (i==res):
            dZdv = np.dot(dZdv,M_i)
        else:
            dZdv = np.dot(dZdv,M[res])
    dZdv = np.dot(dZdv,vecf)
    return (v/Z)*dZdv

def get_Nl(N_res,M,lnv):
    Nl = 0.
    for res in range(0,N_res):
        v = np.exp(lnv[res])
        Nl += calc_Nl_i(N_res,M,res,v)
    return Nl

def get_l(N_res,M,lnv):
    Nl = 0.
    for res in range(1,N_res-1):
        v = np.exp(lnv[res])
        Nl += calc_Nl_i(N_res,M,res,v)
    return Nl


N_res = lnw.shape[2]
M_max = wv_to_M(N_res,lnw[proc_max,samp_max,:][0],lnv[proc_max,samp_max,:][0])
h_i = get_h(N_res,M_max,lnw[proc_max,samp_max,:][0])
h_i_sim = (Nw)/Nk
print h_i
print h_i_sim

mle = pyemma.msm.markov_model(Tmle)
mu = mle.eigenvectors_left(k=1)

Ns = get_Ns(N_res,M_max,lnw[proc_max,samp_max,:][0],lnv[proc_max,samp_max,:][0])
Ns_sim = np.sum(mu*Ns_weight)
Nhb_sim = np.sum(mu*Nhb_weight)
print Ns
print Ns_sim
print Nhb_sim

Nl = get_Nl(N_res,M_max,lnv[proc_max,samp_max,:][0])
Nl_sim = np.sum(mu*Nl_weight)
print Nl
print Nl_sim

l = get_l(N_res,M_max,lnv[proc_max,samp_max,:][0])
l_sim = np.sum(mu*l_weight)
print l
print l_sim

# Now, let's make a kinetic zipper model
k0 = 1.
w_opt = np.exp(lnw[proc_max,samp_max,:][0])[0]
v_opt = np.exp(lnv[proc_max,samp_max,:][0])[0]
hin = 1.*v_opt
kon = k0*hin
K = np.zeros(shape=(len(cc_full),len(cc_full)))
for indi,si in enumerate(cc_full):
    # get the statistical weight for state wi
    LR_ws = state_LR_weights[indi]
    wi = 1.
    for LR_w in LR_ws:
        if ( LR_w == 1 ):
            wi *= v_opt
        elif ( LR_w == 2 ):
            wi *= w_opt
        else:
            wi *= 1.        
    for indj,sj in enumerate(cc_full):
        # only connections between states that differ by a single res state
        if ( np.sum(np.abs(si-sj)) != 1 ):
            continue
        # get the statistical weight for state wj
        LR_ws = state_LR_weights[indj]
        wj = 1.
        for LR_w in LR_ws:
            if ( LR_w == 1 ):
                wj *= v_opt
            elif ( LR_w == 2 ):
                wj *= w_opt
            else:
                wj *= 1.
            
        ind = np.where(np.abs(si-sj)==1)[0][0]
        if ( sj[ind] == 0 ): # flip from c to h
            K[indi,indj] = 1.*kon
        else: # flip from h to c
            K[indi,indj] = (wj/wi)*kon
# now set the diagonal elements
for state in range(len(cc_full)):
    K[state,state] = -np.sum(K[state])

# Get the transition probabilities
T_zip = scipy.linalg.expm(k0*K)
np.sum(T_zip,axis=1)

np.savez('LR_opt',lnw_opt=lnw[proc_max,samp_max,:][0],lnv_opt=lnv[proc_max,samp_max,:][0])
np.savez('eqm_prop',h_i_LR=h_i,h_i_sim=h_i_sim,Ns_LR=Ns,Ns_sim=Ns_sim,Nl_LR=Nl,Nl_sim=Nl_sim,l_LR=l,l_sim=l_sim)
np.savez('kin_zip',K_zip=K,T_zip=T_zip)



