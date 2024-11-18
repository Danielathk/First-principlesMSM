#!/usr/bin/env pypy

from msmbuilder.featurizer import AtomPairsFeaturizer
from msmbuilder.cluster import KCenters
import numpy as np
from msmbuilder.decomposition import tICA
from msmbuilder.dataset import dataset
from matplotlib import pyplot as plt
import matplotlib as mp
import os
from msmbuilder.tpt import net_fluxes
import sys

from msmbuilder.io import load_trajs, save_trajs, save_generic
from msmbuilder.cluster import MiniBatchKMeans
from msmbuilder.cluster import KCenters
from sklearn import preprocessing
from msmbuilder.decomposition import tICA, PCA
from msmbuilder.preprocessing import RobustScaler

from msmbuilder.cluster import KCenters
from matplotlib import pyplot as plt

mp.use('Agg')
k16=16
trajs_16 = []
for i in range(k16):
    temp = np.load('16CO2_ori/outfile-short-a-C%02d.npy'%i)
    trajs_16.append(temp)
    temp = np.load('16CO2_ori/outfile-short-b-C%02d.npy'%i)
    trajs_16.append(temp)
    temp = np.load('16CO2_ori/outfile-short-c-C%02d.npy'%i)
    trajs_16.append(temp)


k=96
trajs = []

#for i in range(k16):
 #   temp = np.load('16CO2_ori/outfile-short-a-C%02d.npy'%i)
 #   trajs.append(temp)
 #   temp = np.load('16CO2_ori/outfile-short-b-C%02d.npy'%i)
 #   trajs.append(temp)
 #   temp = np.load('16CO2_ori/outfile-short-c-C%02d.npy'%i)
 #   trajs.append(temp)

for i in range(k):
    temp = np.load('outfile-short-a-C%02d.npy'%i)
    trajs.append(temp)
    temp = np.load('outfile-short-b-C%02d.npy'%i)
    trajs.append(temp)
    temp = np.load('outfile-short-c-C%02d.npy'%i)
    trajs.append(temp)
#k=k+k16

#O H C
#180/3=60, C:2*4*3=8*3=24;H:8*4*3=32*3=96;O:5*4*3=20*3=60.
#trajs_s = []
#print(trajs)
#print(trajs[0].shape)
#for i in range(3*k):
#    trajs_s.append(trajs[i][::10,:])
#print(trajs[0].shape)
print('performing pca')
pca = PCA(n_components=5)
pca.fit(trajs_16)
pca_trajs = pca.transform(trajs)

np.savetxt('pca_components_.dat',pca.components_)
for i in range(len(pca_trajs)):
    np.savetxt('pca_trajectories_%d'%i, pca_trajs[i])
traj_new=[]
frame0=np.hstack((trajs[0][0],trajs[1][0],trajs[2][0]))


for i in range(k):
    trajsnew=[]
    for j in range(len(trajs[i*3])):
        astack = np.hstack((trajs[i*3][j],trajs[i*3+1][j],trajs[i*3+2][j]))
        trajsnew.append(astack)
    trajsnew=np.array(trajsnew)
    print("i:",trajsnew.shape)
    np.save("outfile-sorted-C%02d"%(i),trajsnew)



txx = np.concatenate(pca_trajs)
clusterer = KCenters(n_clusters=1000,random_state=7)
clustered_trajs = clusterer.fit_transform(pca_trajs)
ktraj_dir = 'ktrajs-extracted-kcenters-%d-%d'%(3,500)

try:
    if not os.path.exists(ktraj_dir):
        os.makedirs(ktraj_dir)
except OSError:
    print("errors while make directory: %s"%ktraj_dir)
    sys.exit()

for i in range(3*k):
    np.save("%s/%08d.npy"%(ktraj_dir,i),clustered_trajs[i],allow_pickle=True)
#save_generic(clusterer, ktraj_pkl)


from msmbuilder.msm import MarkovStateModel, implied_timescales
#data=dataset('cluster',mode='r',fmt='dir-npy',verbose=True)
lag_times=range(200,20000,200)
msm_timescales = implied_timescales(clustered_trajs, lag_times, n_timescales=5,msm=MarkovStateModel(lag_time=250,reversible_type='transpose',ergodic_cutoff='off'))
np.savetxt('msm_timescales_1.txt',msm_timescales)

from msmbuilder.msm import MarkovStateModel
print('building MSM')

msm = MarkovStateModel(lag_time=200, reversible_type='transpose',ergodic_cutoff='off')
msm.fit(clustered_trajs)
print('finished MSM')
from msmbuilder.lumping import PCCAPlus
pcca = PCCAPlus.from_msm(msm, n_macrostates=4)
macro_trajs = pcca.transform(clustered_trajs)
np.savetxt('pcca_2.txt', pcca.microstate_mapping_)
print('finished pcca+')

plt.hexbin(txx[:,0], txx[:,1],bins='log', cmap='viridis')
plt.scatter(clusterer.cluster_centers_[msm.state_labels_, 0],
            clusterer.cluster_centers_[msm.state_labels_, 1],
            s=20,
            c=pcca.microstate_mapping_,
)
plt.xlabel('PC 1')
plt.ylabel('PC 2')
plt.savefig('5macro.png')
#np.savetxt('pcca_2.txt', pcca.microstate_mapping_)



