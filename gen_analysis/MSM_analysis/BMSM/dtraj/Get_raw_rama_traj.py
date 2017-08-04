
# coding: utf-8

# CG MSM from Qhel directly via rama-plot
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
# get the indir from the current dir
cdir = os.getcwd()
flag = False
wdir = 'MSM_analysis'
dir_len = len(wdir)
ctr = 0
while (not flag):
    if ( cdir[ctr:dir_len+ctr] == wdir ):
        flag = True
        dir_ind = ctr
    ctr += 1

indir = cdir[:ctr+dir_len-(len(wdir)+1)]

#indir = '/data/isilon/rudzinski/cluster_tmp/ALA7/ALA7_hybrid_AMBER_Go/wDB-HP_inter/Find_Tfold/NC_CA/2016_10_11/epsNC-8/epsNC-8_epsdb-0.00epsNC_epshp-0.3epsNC/T-550/'

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

dtraj_rama = []
dtraj_rama.append( dtraj_rama_2 )
dtraj_rama.append( dtraj_rama_3 )
dtraj_rama.append( dtraj_rama_4 )
dtraj_rama.append( dtraj_rama_5 )
dtraj_rama.append( dtraj_rama_6 )



np.save('traj_rama.npy', dtraj_rama)
