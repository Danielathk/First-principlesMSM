#!/usr/bin/env pypy

'''
microstate clustering, here we need to assemble 3 characters of oxygen to represent the state of central carbon atom
'''


from msmbuilder.featurizer import AtomPairsFeaturizer
import numpy as np
from msmbuilder.decomposition import tICA
from msmbuilder.dataset import dataset
from matplotlib import pyplot as plt
import matplotlib as mp
import os
from msmbuilder.tpt import net_fluxes
import sys

from msmbuilder.io import load_trajs, save_trajs, save_generic
from sklearn.cluster import MiniBatchKMeans, KMeans
from msmbuilder.cluster import KCenters
from sklearn import preprocessing
from msmbuilder.decomposition import tICA, PCA


mp.use('Agg')

ntrajs=96
dim=2
nclusters=90
ktraj_dir = 'ktrajs-extracted-kmeans-%d-%d'%(dim, nclusters)

clustered_trajs = []
for i in range(3*ntrajs):
    temp = np.load('%s/%06d.npy'%(ktraj_dir,i))
    clustered_trajs.append(temp)

clustered_trajs = np.array(clustered_trajs)
#txx= np.concatenate(clustered_trajs)

#np.save('cluster_center',clusters.cluster_centers_)
#np.save('cluster_label',clusters.labels_)

#clustered_trajs =[]
#for i in range(3*ntrajs):
#    clustered_trajs.append(clusters.predict(pcatrajs[i]))
#clustered_trajs=np.array(clustered_trajs)

#ktraj_dir = 'ktrajs-extracted-kmeans-%d-%d'%(dim, nclusters)

#try: 
#    if not os.path.exists(ktraj_dir):
#        os.makedirs(ktraj_dir)
#except OSError:
#    print("errors while make directory: %s"%ktraj_dir)

#for i in range(3*ntrajs):
#    np.save("%s/%06d.npy"%(ktraj_dir,i),clustered_trajs[i],allow_pickle=True)

cluster_center = np.load('cluster_center.npy')
mapping_a=[]
for i in range(nclusters):
    if cluster_center[i][1]<0.2:
        if cluster_center[i][0]<0.1:
            mapping_a.append(0)
        if cluster_center[i][0]>=0.1:
            mapping_a.append(1)
    if cluster_center[i][1]>=0.2:
        mapping_a.append(2)
print(mapping_a)
map = np.empty(nclusters)
map[:] = np.NaN
rec = 2
for i in range(nclusters):
    if mapping_a[i] !=2:
        map[i]=mapping_a[i]
    else:
        map[i]=rec
        rec=rec+1
print(map, rec)
clustered_new=[]
for i in range(3*ntrajs):
    temp = np.zeros(len(clustered_trajs[i]))
    for j in range(len(clustered_trajs[i])):
        temp[j]=map[clustered_trajs[i][j]]
    clustered_new.append(temp)
clustered_new = np.array(clustered_new)

ktrajnew_dir = 'ktrajs-new-kmeans-%d-%d'%(dim, rec)
try:
    if not os.path.exists(ktrajnew_dir):
        os.makedirs(ktrajnew_dir)
except OSError:
    print("errors while make directory: %s"%ktrajnew_dir)

for i in range(3*ntrajs):
    np.save("%s/%06d.npy"%(ktrajnew_dir,i),clustered_new[i],allow_pickle=True)


txx= np.concatenate(clustered_new)

micro_trajs = []
for i in range(ntrajs):
    mt=[]
    for j in range(len(clustered_new[3*i])):
        ms = []
        ms.append(int(clustered_new[3*i][j]))
        ms.append(int(clustered_new[3*i+1][j]))
        ms.append(int(clustered_new[3*i+2][j]))
        ms.sort()
        mt.append(ms[0]*rec*rec+ms[1]*rec+ms[2])
    mt=np.array(mt)
    micro_trajs.append(mt)
micro_trajs=np.array(micro_trajs)

m_trajs_C = np.concatenate(micro_trajs)
m_trajs_C.sort()

#print(max(m_trajs_C))
mapping =[]
mapping_back = np.empty(int(max(m_trajs_C)+1))
mapping_back[:] = np.NaN
mapping.append(m_trajs_C[0])
mapping_back[int(m_trajs_C[0])] = len(mapping)-1
for i in range(len(m_trajs_C)):
    if (m_trajs_C[i]> mapping[len(mapping)-1]):
        mapping.append(m_trajs_C[i])
        mapping_back[int(m_trajs_C[i])] = len(mapping)-1

np.save('mapping.dat',mapping)
np.save('map_back.dat',mapping_back)
print(max(mapping_back))

micro_trajs_new = []
for i in range(ntrajs):
    temp = []
    for j in range(len(micro_trajs[i])):
        temp.append(mapping_back[micro_trajs[i][j]])
    temp = np.array(temp)
    micro_trajs_new.append(temp)
micro_trajs_new = np.array(micro_trajs_new)



microtraj_dir = 'micro_trajs-kmeans-%d-%d'%(dim, rec)

try:
    if not os.path.exists(microtraj_dir):
        os.makedirs(microtraj_dir)
except OSError:
    print("errors while make directory: %s"%microtraj_dir)

for i in range(ntrajs):
    np.save("%s/%06d.npy"%(microtraj_dir,i),micro_trajs_new[i],allow_pickle=True)



