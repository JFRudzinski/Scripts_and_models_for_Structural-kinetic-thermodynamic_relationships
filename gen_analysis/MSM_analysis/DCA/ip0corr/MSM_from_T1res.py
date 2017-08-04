# imports
import numpy as np
from sympy import *
import pyemma
pyemma.__version__
import os
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
import msmbuilder
from msmbuilder.msm.ratematrix import ContinuousTimeMSM
import scipy
from msmtools.analysis.dense.decomposition import eigenvectors, eigenvalues
import operator
from scipy.optimize import fsolve, minimize
from copy import deepcopy

# for saving objects
import pickle
def save_object(filename, obj):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

# user input
Nres = 5
Ncorr = 0 # == k
T_basenm = 'T1res'
flag_seqdep = False
flag_bndry_cond = False
bndry_cond_path = ''
tol = 1e-2
max_steps = 10
end_state = 1 # 1 betas / 0 helices at the ends

# fixed variables
S2 = [0,1]

# Read in system-specific quantities
Pkp1res_inp = []
Nsets = Nres-Ncorr
if ( flag_seqdep ):
    Pkp1res_inp = np.load(T_basenm+'_seqdep.npy')
else:
    for group in range(Nsets):
        Pkp1res_inp.append(np.load(T_basenm+'.npy'))
cckp1res = np.load('../cc1res.npy')
cckres = np.array([[S2[0]],[S2[1]]]) #np.array( S2 )  # np.load('cc3res.npy')
cc_full = np.load('../cc_full.npy')

# Some useful functions
def add_Pkp1res_stat_dict_vals( dic, Pkp1res_stat, Pkp1res_stat_val, states_kp1res ):
    for res in range(Nres-Ncorr):
        for state in range(len(states_kp1res)):
            var = Pkp1res_stat[res][state]
            dic[var] = Pkp1res_stat_val[res][state]
    return dic

def sub_inp( eqns, dic ):
    for eqn in range(len(eqns)):
        eqns[eqn] = eqns[eqn].subs(dic)
    return eqns

def lambdify_vec( eqns, var ):
    fncs = [lambdify((var), eqn, modules='numpy') for eqn in eqns]
    return fncs

def gen_fv( fncs ):
    return lambda zz: np.array([fnc(*zz) for fnc in fncs])

def gen_jv( fncs ):
    return lambda zz: np.array([ [ifnc(*zz) for ifnc in jfncs] for jfncs in fncs ])

def init_soln( N ):
    return 0.5*np.ones(N)

def sumsq_eqns( eqns, var ):
    sumsq = 0.
    for eqn in range(len(eqns)):
        sumsq += eqns[eqn].subs(var)**2
    return sumsq

# Set up the equations for the kres cond probs
def get_poss_states( Nres ):
    states = []
    state = [0]*Nres
    for gstate in range(2**Nres):
        states.append(state[:])
        state[0] = (state[0]+1) % 2
        for res in range(1,Nres):
            if ( (state[res-1]+1) % 2 == 1 ):
                state[res] = (state[res]+1) % 2
            else:
                break
    return states

def get_poss_states_endignore( Nres, Nend_front, Nend_back, end_state ):
    states = []
    state = [0]*Nres
    for gstate in range(2**Nres):
        state[0] = (state[0]+1) % 2
        for res in range(1,Nres):
            if ( (state[res-1]+1) % 2 == 1 ):
                state[res] = (state[res]+1) % 2
            else:
                break

        flag_add = True
        for res in range(Nend_front):
            if ( state[res] != end_state ):
                flag_add = False
        for res in range(Nend_back):
            if ( state[-(res+1)] != end_state ):
                flag_add = False

        if ( flag_add ):
            states.append(state[:])
     
    return states

# Directly calculate the approximate TN matrix
# first, get all the possible states of the full system
states_Nres = cc_full.astype(int)
TN = np.ones(shape=(len(states_Nres),len(states_Nres)))
for si in range(len(states_Nres)):
    for sj in range(len(states_Nres)):
        for res in range(Nres):
            inpi = states_Nres[si][res]
            inpj = states_Nres[sj][res]
            TN[si][sj] *= Pkp1res_inp[res][inpi][inpj]
        

np.save('TN_opt',TN)




