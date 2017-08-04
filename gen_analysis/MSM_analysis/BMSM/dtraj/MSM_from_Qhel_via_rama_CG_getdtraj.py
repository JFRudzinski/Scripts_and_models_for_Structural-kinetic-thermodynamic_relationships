
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

# In[4]:

indir = '/data/isilon/rudzinski/cluster_tmp/ALA7/ALA7_hybrid_AMBER_Go/wDB-HP_inter/Find_Tfold/NC_CA/2016_10_11/epsNC-11/epsNC-11_epsdb-0.05epsNC_epshp-0.6epsNC/T-250/'

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

for i in range( len(dtraj_rama_2) ):
    dtraj_rama_2[i][np.where(dtraj_rama_2[i][:,1] < -125)[0],1] += 360  
    dtraj_rama_3[i][np.where(dtraj_rama_3[i][:,1] < -125)[0],1] += 360
    dtraj_rama_4[i][np.where(dtraj_rama_4[i][:,1] < -125)[0],1] += 360
    dtraj_rama_5[i][np.where(dtraj_rama_5[i][:,1] < -125)[0],1] += 360
    dtraj_rama_6[i][np.where(dtraj_rama_6[i][:,1] < -125)[0],1] += 360


# In[6]:

dtraj_phi_2 = []
dtraj_phi_3 = []
dtraj_phi_4 = []
dtraj_phi_5 = []
dtraj_phi_6 = []
for i in range( len(dtraj_rama_2) ):
    dtraj_phi_2.append(dtraj_rama_2[i][:,1])
    dtraj_phi_3.append(dtraj_rama_3[i][:,1])
    dtraj_phi_4.append(dtraj_rama_4[i][:,1])
    dtraj_phi_5.append(dtraj_rama_5[i][:,1])
    dtraj_phi_6.append(dtraj_rama_6[i][:,1])


# **simple clustering along psi only for discretization**

# In[7]:

n_clusters = 2     # number of k-means clusters


# In[8]:

clustering_rama_2 = coor.cluster_kmeans(dtraj_phi_2,k=n_clusters,max_iter=100, tolerance=1e-12, fixed_seed=True)
clustering_rama_3 = coor.cluster_kmeans(dtraj_phi_3,k=n_clusters,max_iter=100, tolerance=1e-12, fixed_seed=True)
clustering_rama_4 = coor.cluster_kmeans(dtraj_phi_4,k=n_clusters,max_iter=100, tolerance=1e-12, fixed_seed=True)
clustering_rama_5 = coor.cluster_kmeans(dtraj_phi_5,k=n_clusters,max_iter=100, tolerance=1e-12, fixed_seed=True)
clustering_rama_6 = coor.cluster_kmeans(dtraj_phi_6,k=n_clusters,max_iter=100, tolerance=1e-12, fixed_seed=True)


# In[9]:

cc_rama_2 = clustering_rama_2.clustercenters[:,0]
cc_rama_3 = clustering_rama_3.clustercenters[:,0]
cc_rama_4 = clustering_rama_4.clustercenters[:,0]
cc_rama_5 = clustering_rama_5.clustercenters[:,0]
cc_rama_6 = clustering_rama_6.clustercenters[:,0]
print cc_rama_2
print cc_rama_3
print cc_rama_4
print cc_rama_5
print cc_rama_6


# In[10]:

dtrajs_rama_2 = clustering_rama_2.dtrajs
dtrajs_rama_3 = clustering_rama_3.dtrajs
dtrajs_rama_4 = clustering_rama_4.dtrajs
dtrajs_rama_5 = clustering_rama_5.dtrajs
dtrajs_rama_6 = clustering_rama_6.dtrajs


# In[11]:

for traj in range( len(dtraj_rama_2) ):
    if ( dtrajs_rama_2[traj][np.where(dtraj_phi_2[traj] < 0)[0][0]] != 0 ):
        dtrajs_rama_2[traj][np.where(dtrajs_rama_2[traj] == 0)[0]] -= 1
        dtrajs_rama_2[traj][np.where(dtrajs_rama_2[traj] == 1)[0]] -= 1
        dtrajs_rama_2[traj][np.where(dtrajs_rama_2[traj] == -1)[0]] += 2
    if ( dtrajs_rama_3[traj][np.where(dtraj_phi_3[traj] < 0)[0][0]] != 0 ):
        dtrajs_rama_3[traj][np.where(dtrajs_rama_3[traj] == 0)[0]] -= 1
        dtrajs_rama_3[traj][np.where(dtrajs_rama_3[traj] == 1)[0]] -= 1
        dtrajs_rama_3[traj][np.where(dtrajs_rama_3[traj] == -1)[0]] += 2
    if ( dtrajs_rama_4[traj][np.where(dtraj_phi_4[traj] < 0)[0][0]] != 0 ):
        dtrajs_rama_4[traj][np.where(dtrajs_rama_4[traj] == 0)[0]] -= 1
        dtrajs_rama_4[traj][np.where(dtrajs_rama_4[traj] == 1)[0]] -= 1
        dtrajs_rama_4[traj][np.where(dtrajs_rama_4[traj] == -1)[0]] += 2
    if ( dtrajs_rama_5[traj][np.where(dtraj_phi_5[traj] < 0)[0][0]] != 0 ):
        dtrajs_rama_5[traj][np.where(dtrajs_rama_5[traj] == 0)[0]] -= 1
        dtrajs_rama_5[traj][np.where(dtrajs_rama_5[traj] == 1)[0]] -= 1
        dtrajs_rama_5[traj][np.where(dtrajs_rama_5[traj] == -1)[0]] += 2
    if ( dtrajs_rama_6[traj][np.where(dtraj_phi_6[traj] < 0)[0][0]] != 0 ):
        dtrajs_rama_6[traj][np.where(dtrajs_rama_6[traj] == 0)[0]] -= 1
        dtrajs_rama_6[traj][np.where(dtrajs_rama_6[traj] == 1)[0]] -= 1
        dtrajs_rama_6[traj][np.where(dtrajs_rama_6[traj] == -1)[0]] += 2


# In[12]:

dtrajs = []
for i in range( len(dtraj_rama_2) ):
    dtrajs.append( np.vstack( (dtrajs_rama_2[i], dtrajs_rama_3[i], dtrajs_rama_4[i], dtrajs_rama_5[i], dtrajs_rama_6[i]) ).T )
    dtrajs[i].astype('int64')


# In[18]:

dtrajs_sum = []
for i in range( len(dtrajs) ):
    dtrajs_sum.append( np.sum(dtrajs[i],axis=1) )
    dtrajs_sum[i].astype('int64')


# In[19]:

# we need a single dimensional identifier of the microstate, can we cluster to automize?
n_clusters = 6
clustering = coor.cluster_regspace(dtrajs_sum,max_centers=n_clusters,dmin=0.5)
#clustering = coor.cluster_kmeans(dtrajs,k=n_clusters,max_iter=100, tolerance=1e-12, fixed_seed=True)


# In[20]:

dtrajs_1D = clustering.dtrajs


# In[21]:

cc = clustering.clustercenters[:]
cc


# In[22]:

sorted_list = sorted(cc)
cc_sorted = np.array(sorted_list)[:]
cc_sorted = cc_sorted.tolist()
for i in range(len(cc)):
    cc_sorted[i] = map(int,cc_sorted[i])
cc_sorted


# In[23]:

for traj in range(len(dtrajs_1D)):
    dtrajs_1D[traj] = cc[dtrajs_1D[traj]].astype(int)


# In[24]:

# we need a single dimensional identifier of the microstate, can we cluster to automize?
n_clusters = 32
clustering = coor.cluster_regspace(dtrajs,max_centers=n_clusters,dmin=0.5)
#clustering = coor.cluster_kmeans(dtrajs,k=n_clusters,max_iter=100, tolerance=1e-12, fixed_seed=True)


# In[25]:

dtrajs_full_1D = clustering.dtrajs


# In[26]:

cc = clustering.clustercenters[:]
cc


# In[27]:

# calculate the number of adjacent and almost adjacent alphas for each state
cc_nn = np.zeros(len(cc))
cc_nnn = np.zeros(len(cc))
for i in range(len(cc)):
    if ( np.where(cc[i]==1)[0].shape[0] == 2):
        d1 = np.where(cc[i]==1)[0][1] - np.where(cc[i]==1)[0][0]
        if (d1==1):
            cc_nn[i] += 1
        elif (d1==2):
            cc_nnn[i] += 1
    elif ( np.where(cc[i]==1)[0].shape[0] == 3):
        d1 = np.where(cc[i]==1)[0][1] - np.where(cc[i]==1)[0][0]
        d2 = np.where(cc[i]==1)[0][2] - np.where(cc[i]==1)[0][1]
        if (d1==1):
            cc_nn[i] += 1
        elif (d1==2):
            cc_nnn[i] += 1
        if (d2==1):
            cc_nn[i] += 1
        elif (d2==2):
            cc_nnn[i] += 1
    elif ( np.where(cc[i]==1)[0].shape[0] == 4):
        d1 = np.where(cc[i]==1)[0][1] - np.where(cc[i]==1)[0][0]
        d2 = np.where(cc[i]==1)[0][2] - np.where(cc[i]==1)[0][1]
        d3 = np.where(cc[i]==1)[0][3] - np.where(cc[i]==1)[0][2]
        if (d1==1):
            cc_nn[i] += 1
        elif (d1==2):
            cc_nnn[i] += 1
        if (d2==1):
            cc_nn[i] += 1
        elif (d2==2):
            cc_nnn[i] += 1
        if (d3==1):
            cc_nn[i] += 1
        elif (d3==2):
            cc_nnn[i] += 1


# In[28]:

lcc = np.arange(len(cc))


# In[29]:

cc_stack = []
for i in range(len(cc)):
    cc_stack.append(np.hstack((np.sum(cc,axis=1)[i],cc_nn[i],cc_nnn[i],cc[i],lcc[i])))


# In[30]:

sorted_list = sorted(cc_stack, key=operator.itemgetter(0,1,2,5,4,6,3))


# In[31]:

cc_sorted = np.array(sorted_list)[:,3:8]
cc_sorted = cc_sorted.tolist()
for i in range(len(cc)):
    cc_sorted[i] = map(int,cc_sorted[i])
cc_sorted


# In[32]:

lcc_sorted = np.array(sorted_list)[:,8]
lcc_sorted = lcc_sorted.tolist()
lcc_sorted = map(int, lcc_sorted)
lcc_sorted


# In[33]:

lcc_sorted_map = []
for i in range(cc.shape[0]):
    lcc_sorted_map.append( np.where(lcc_sorted==lcc[i])[0][0] )
lcc_sorted_map = np.array(lcc_sorted_map)
print lcc_sorted_map


# In[34]:

dtrajs_full_1D


# In[35]:

dtrajs_full_1D_new = []
for traj in range(len(dtrajs_full_1D)):
    dtrajs_full_1D_new.append( lcc_sorted_map[dtrajs_full_1D[traj]].astype(int) )


# In[36]:

np.save('dtrajs_rama_2st_res2', dtrajs_rama_2)
np.save('dtrajs_rama_2st_res3', dtrajs_rama_3)
np.save('dtrajs_rama_2st_res4', dtrajs_rama_4)
np.save('dtrajs_rama_2st_res5', dtrajs_rama_5)
np.save('dtrajs_rama_2st_res6', dtrajs_rama_6)
#
np.save('dtrajs_rama_2st_allres', dtrajs)
#
np.save('dtrajs_rama_2st_allres_1D', dtrajs_full_1D_new)
#
np.save('dtrajs_Qhel_allres', dtrajs_sum)
#
np.save('dtrajs_Qhel_allres_1D', dtrajs_1D)


# In[ ]:




# In[37]:

tau = 25
nts = n_clusters
Cmat = pyemma.msm.estimation.count_matrix(dtrajs_full_1D_new, tau, sliding=True, sparse_return=False)
lcc = pyemma.msm.estimation.largest_connected_set(Cmat, directed=True)
Cmat = pyemma.msm.estimation.largest_connected_submatrix(Cmat, directed=True, lcc=lcc)
Tmle = pyemma.msm.estimation.transition_matrix(Cmat, reversible=True)


# In[38]:

mle = pyemma.msm.markov_model(Tmle)
evals_mle = mle.eigenvalues(k=n_clusters)
evecs_mle = mle.eigenvectors_left(k=n_clusters)


# In[39]:

