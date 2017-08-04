
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


# Read in the data
# ------

# The LR sampling output data
Nproc = 16
lnw = []
lnv = []
lnc = []
lnn = []
L = []
for i in range(Nproc):
    data = np.load('output_variables_'+str(i)+'.npz')
    lnw.append(data['lnw'])
    lnv.append(data['lnv'])
    lnc.append(data['lnc'])
    lnn.append(data['lnn'])
    L.append(data['L'])
lnw = np.array(lnw)
lnv = np.array(lnv)
lnc = np.array(lnc)
lnn = np.array(lnn)
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
#
LR_const_params_cap = np.load('../LR_const_params_cap.npz')
Nv_cap = LR_const_params_cap['Nv']
Nw_cap = LR_const_params_cap['Nw']
Nc_cap = LR_const_params_cap['Nc']
Nn_cap = LR_const_params_cap['Nn']

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
state_LR_weights_cap = np.load(indir+'state_LR_weights_cap.npy')

# The MSM
Tmle = np.load('../../MLE/T_mle.npy')

proc_max,samp_max = np.where(L==np.max(L))
proc_avg,samp_avg = np.where( L-np.mean(L)==np.min(np.abs(L-np.mean(L))) )


def wv_to_M(N_res,lnw,lnv,lnc,lnn):
    M = np.zeros(shape=(N_res,6,6))
    for res in range(N_res):
        w = np.exp(lnw[res])
        v = np.exp(lnv[res])
        if ( res+2 < N_res ):
            c = np.exp(lnc[res+2])
        else:
            c = 0.
        if ( res-2 >= 0 ):
            n = np.exp(lnn[res-2])
        else:
            n = 0.        
        M[res,0,0] = w
        M[res,0,1] = w*c
        M[res,1,3] = v
        M[res,2,0] = w*n
        M[res,2,1] = w*n*c
        M[res,3,4] = 1.
        M[res,3,5] = 1.
        M[res,4,2] = v
        M[res,4,3] = v
        M[res,5,4] = 1.
        M[res,5,5] = 1.
        
    return M

def M_to_Z(N_res,M):
    vec0 = np.array([0.,0.,0.,0.,0.,1.])
    vecf = np.array([0.,0.,0.,0.,0.,1.])
    Z = np.dot(vec0,M[0])
    for res in range(1,N_res):
        Z = np.dot(Z,M[res])
    Z = np.dot(Z,vecf)
    return Z

def calc_h_i(N_res,M,i,w):
    from copy import deepcopy
    M_i = np.zeros(shape=M[i].shape)
    M_i[0,0] = deepcopy(M[i][0,0]) / w
    M_i[0,1] = deepcopy(M[i][0,1]) / w
    M_i[2,0] = deepcopy(M[i][2,0]) / w
    M_i[2,1] = deepcopy(M[i][2,1]) / w
      
    vec0 = np.array([0.,0.,0.,0.,0.,1.])
    vecf = np.array([0.,0.,0.,0.,0.,1.])
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

def calc_Ns_i(N_res,M,i,w,c):
    # we are calculating d^2 lnZ / (dlnw_i-1 dlnvi)
    from copy import deepcopy
    dMidw = np.zeros(shape=M[i].shape)
    dMidw[0,0] = deepcopy(M[i][0,0]) / w
    dMidw[0,1] = deepcopy(M[i][0,1]) / w
    dMidw[2,0] = deepcopy(M[i][2,0]) / w
    dMidw[2,1] = deepcopy(M[i][2,1]) / w
    dMidc = np.zeros(shape=M[i].shape)
    dMidc[0,1] = deepcopy(M[i][0,1]) / c
    dMidc[2,1] = deepcopy(M[i][2,1]) / c
    d2Mi = np.zeros(shape=M[i].shape)
    d2Mi[0,1] = deepcopy(M[i][0,1]) / (w*c)
    d2Mi[2,1] = deepcopy(M[i][2,1]) / (w*c)    
    
    vec0 = np.array([0.,0.,0.,0.,0.,1.])
    vecf = np.array([0.,0.,0.,0.,0.,1.])
    Z = M_to_Z(N_res,M)
    # dZ/dlnw
    if ( i == 0 ):
        dZdlnw = np.dot(vec0,dMidw)    
    else:
        dZdlnw = np.dot(vec0,M[0])
    for res in range(1,N_res):
        if (i==res):
            dZdlnw = np.dot(dZdlnw,dMidw)
        else:
            dZdlnw = np.dot(dZdlnw,M[res])
    dZdlnw = np.dot(dZdlnw,vecf)
    dZdlnw = w*dZdlnw
    # dZ/dlnc
    if ( i == 0 ):
        dZdlnc = np.dot(vec0,dMidc)    
    else:
        dZdlnc = np.dot(vec0,M[0])
    for res in range(1,N_res):
        if (i==res):
            dZdlnc = np.dot(dZdlnc,dMidc)
        else:
            dZdlnc = np.dot(dZdlnc,M[res])
    dZdlnc = np.dot(dZdlnc,vecf)
    dZdlnc = c*dZdlnc
    # d2Z/dlnwlnc
    if ( i == 0 ):
        d2Z = np.dot(vec0,d2Mi)    
    else:
        d2Z = np.dot(vec0,M[0])
    for res in range(1,N_res):
        if (i==res):
            d2Z = np.dot(d2Z,d2Mi)
        else:
            d2Z = np.dot(d2Z,M[res])
    d2Z = np.dot(d2Z,vecf)
    d2Z = w*c*d2Z
    
    return (1./(Z**2))*( (Z*d2Z) - (dZdlnw*dZdlnc) )

def get_Ns(N_res,M,lnw,lnc):
    Ns = 0.
    for res in range(2,N_res-2):
        w = np.exp(lnw[res])
        c = np.exp(lnc[res+2])
        Ns += calc_Ns_i(N_res,M,res,w,c)
    return Ns

def calc_Nl_i(N_res,M,i,v):
    from copy import deepcopy
    M_i = np.zeros(shape=M[i].shape)
    M_i[4,3] = deepcopy(M[i][4,3])/v
    vec0 = np.array([0.,0.,0.,0.,0.,1.])
    vecf = np.array([0.,0.,0.,0.,0.,1.])
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
    for res in range(1,N_res-1):
        v = np.exp(lnv[res])
        Nl += calc_Nl_i(N_res,M,res,v)
    return Nl

def get_l(N_res,M,lnv):
    Nl = 0.
    for res in range(2,N_res-2):
        v = np.exp(lnv[res])
        Nl += calc_Nl_i(N_res,M,res,v)
    return Nl


N_res = lnw.shape[2]
M_max = wv_to_M(N_res,lnw[proc_max,samp_max,:][0],lnv[proc_max,samp_max,:][0],lnc[proc_max,samp_max,:][0],lnn[proc_max,samp_max,:][0])
h_i = get_h(N_res,M_max,lnw[proc_max,samp_max,:][0])
h_i_sim = (Nw_cap)/Nk
print h_i
print h_i_sim

mle = pyemma.msm.markov_model(Tmle)
mu = mle.eigenvectors_left(k=1)

Ns = get_Ns(N_res,M_max,lnw[proc_max,samp_max,:][0],lnc[proc_max,samp_max,:][0])
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
c_opt = np.exp(lnc[proc_max,samp_max,:][0])[N_res-1]
n_opt = np.exp(lnn[proc_max,samp_max,:][0])[0]
hin = 1.*v_opt
kon = k0*hin
K = np.zeros(shape=(len(cc_full),len(cc_full)))
for indi,si in enumerate(cc_full):
    # get the statistical weight for state wi
    LR_ws = state_LR_weights_cap[indi,1:N_res-1]
    wi = 1.
    for LR_ind,LR_w in enumerate(LR_ws):
        # first get the capped states
        if ( LR_ind+2 < N_res ):
            cip2 = state_LR_weights_cap[indi,LR_ind+2]
        else:
            cip2 = 0
        if ( LR_ind-2 >= 0 ):
            nim2 = state_LR_weights_cap[indi,LR_ind-2]
        else:
            nim2 = 0
        # now for the normal weights
        if ( LR_w == 1 ):
            wi *= v_opt
        elif ( LR_w == 2 ):
            wi *= w_opt
        #else:
        #    wi *= 1.
        # and the additional capping effect
        if ( cip2 == 4 ):
            if ( LR_ind+2 == N_res-1 ):
                wi *= c_opt
        if ( nim2 == 3 ):
            if ( LR_ind-2 == 0 ):
                wi *= n_opt
            
    for indj,sj in enumerate(cc_full):
        # only connections between states that differ by a single res state
        if ( np.sum(np.abs(si-sj)) != 1 ):
            continue
        ind = np.where(np.abs(si-sj)==1)[0][0]
        
        # get the statistical weight for state wj
        LR_ws = state_LR_weights_cap[indj,1:N_res-1]
        wj = 1.
        for LR_ind,LR_w in enumerate(LR_ws):
            # first get the capped states
            if ( LR_ind+2 < N_res ):
                cjp2 = state_LR_weights_cap[indj,LR_ind+2]
            else:
                cjp2 = 0
            if ( LR_ind-2 >= 0 ):
                njm2 = state_LR_weights_cap[indj,LR_ind-2]
            else:
                njm2 = 0
            # now for the normal weights
            if ( LR_w == 1 ):
                wj *= v_opt
            elif ( LR_w == 2 ):
                wj *= w_opt
            #else:
            #    wj *= 1.
            # and the additional capping effect
            if ( cjp2 == 4 ):
                if ( LR_ind+2 == N_res-1 ):
                    wj *= c_opt
            if ( njm2 == 3 ):
                if ( LR_ind-2 == 0 ):
                    wj *= n_opt   
        
        # get the rates
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

np.savez('LR_opt',lnw_opt=lnw[proc_max,samp_max,:][0],lnv_opt=lnv[proc_max,samp_max,:][0],lnc_opt=lnc[proc_max,samp_max,:][0],lnn_opt=lnn[proc_max,samp_max,:][0])
np.savez('eqm_prop',h_i_LR=h_i,h_i_sim=h_i_sim,Ns_LR=Ns,Ns_sim=Ns_sim,Nl_LR=Nl,Nl_sim=Nl_sim,l_LR=l,l_sim=l_sim)
np.savez('kin_zip',K_zip=K,T_zip=T_zip)

