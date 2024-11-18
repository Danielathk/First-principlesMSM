#!/usr/bin/env pypy

'''
combine pca of three oxygen atoms that are cloest to center carbon atoms
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
ntrajs_a = 32
ntrajs_b = 48
ntrajs_ori = 16
ntrajs = ntrajs_ori+ ntrajs_a + ntrajs_b
dim=2
nclusters=90

pcatrajs_ori = []
for i in range(3*ntrajs_ori):
    temp = np.loadtxt('../../../16CO2_ori/pca_trajectories_%d'%i)
    pcatrajs_ori.append(temp[:,0:dim])

pcatrajs_ori = np.array(pcatrajs_ori)
#print(pcatrajs.shape)

txx_o= np.concatenate(pcatrajs_ori)

clusters = MiniBatchKMeans(n_clusters = nclusters,random_state=982)
#clusters = KMeans(n_clusters = nclusters,random_state=7)

clusters.fit(txx_o)
np.save('cluster_center',clusters.cluster_centers_)
np.save('cluster_label',clusters.labels_)

pcatrajs = []
for i in range(3*ntrajs_ori):
    temp = np.loadtxt('../../../16CO2_ori/pca_trajectories_%d'%i)
    pcatrajs.append(temp[:,0:dim])

for i in range(3*ntrajs_a):
    temp = np.loadtxt('../../../a/pca_trajectories_%d'%i)
    pcatrajs.append(temp[:,0:dim])
for i in range(3*ntrajs_b):
    temp = np.loadtxt('../../../b/pca_trajectories_%d'%i)
    pcatrajs.append(temp[:,0:dim])
pcatrajs = np.array(pcatrajs)


clustered_trajs =[]
for i in range(3*ntrajs):
    clustered_trajs.append(clusters.predict(pcatrajs[i]))
clustered_trajs=np.array(clustered_trajs)

ktraj_dir = 'ktrajs-extracted-kmeans-%d-%d'%(dim, nclusters)

try: 
    if not os.path.exists(ktraj_dir):
        os.makedirs(ktraj_dir)
except OSError:
    print("errors while make directory: %s"%ktraj_dir)

for i in range(3*ntrajs):
    np.save("%s/%06d.npy"%(ktraj_dir,i),clustered_trajs[i],allow_pickle=True)
txx = np.concatenate(pcatrajs)
reduced_data = txx

np.savetxt('txx_plot.txt', txx)
# Step size of the mesh. Decrease to increase the quality of the VQ.
h = 0.01  # point in the mesh [x_min, x_max]x[y_min, y_max].

# Plot the decision boundary. For that, we will assign a color to each
x_min, x_max = reduced_data[:, 0].min() - 0.01, reduced_data[:, 0].max() + 0.01
y_min, y_max = reduced_data[:, 1].min() - 0.01, reduced_data[:, 1].max() + 0.01
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

# Obtain labels for each point in mesh. Use last trained model.
Z = clusters.predict(np.c_[xx.ravel(), yy.ravel()])

# Put the result into a color plot
Z = Z.reshape(xx.shape)
plt.figure(1)
plt.clf()
plt.imshow(
    Z,
    interpolation="nearest",
    extent=(xx.min(), xx.max(), yy.min(), yy.max()),
    cmap=plt.cm.Paired,
    aspect="auto",
    origin="lower",
)

#plt.plot(reduced_data[:, 0], reduced_data[:, 1], "k.", markersize=2)
# Plot the centroids as a white X
centroids = clusters.cluster_centers_
plt.scatter(
    centroids[:, 0],
    centroids[:, 1],
    marker="x",
    s=169,
    linewidths=3,
    color="w",
    zorder=10,
)
plt.title(
    "K-means clustering on the digits dataset (PCA-reduced data)\n"
    "Centroids are marked with white cross"
)
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.xticks(())
plt.yticks(())
plt.savefig('3.png',dpi=1000)
plt.clf()

#fig=plt.figure()
#cmap = mp.cm.coolwarm
#cmap = 'Purples'
cmap_reversed = mp.cm.get_cmap('RdYlBu_r')
#cmap = 'RdBu'

plt.hexbin(txx[:,0], txx[:,1],bins='log', cmap=cmap_reversed, vmin=0.01)
plt.colorbar()
plt.scatter(clusters.cluster_centers_[:, 0],
            clusters.cluster_centers_[:, 1],
            s=20,
            
)
#cmap = mpl.cm.cool
#norm = mpl.colors.Normalize(vmin=5, vmax=10)

#fig.colorbar()

plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.xticks(([-1, -0.5, 0, 0.5, 1, 1.5, 2]))
plt.yticks(([-0.5, 0, 0.5 , 1, 1.5]))

plt.xlabel('PC 1')
plt.ylabel('PC 2')
plt.savefig('1.png', dpi=700)
plt.clf()




cmap_reversed = mp.cm.get_cmap('RdYlBu_r')
#cmap = 'RdBu'

plt.hexbin(txx[:,0], txx[:,1],bins='log', cmap=cmap_reversed, vmin=0.1)
plt.colorbar()
plt.scatter(clusters.cluster_centers_[:, 0],
            clusters.cluster_centers_[:, 1],
            s=20,
            
)
#cmap = mpl.cm.cool
#norm = mpl.colors.Normalize(vmin=5, vmax=10)

#fig.colorbar()

plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.xticks(([-1, -0.5, 0, 0.5, 1, 1.5, 2]))
plt.yticks(([-0.5, 0, 0.5 , 1, 1.5]))

plt.xlabel('PC 1')
plt.ylabel('PC 2')
plt.savefig('1_b.png', dpi=700)
plt.clf()


plt.hexbin(txx[:,0], txx[:,1],bins='log', cmap=cmap_reversed,vmin=0.1)
plt.xlabel('PC 1')
plt.ylabel('PC 2')
plt.savefig('2.png', dpi=1000)



