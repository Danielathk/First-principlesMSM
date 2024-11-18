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



fr=open('data.lst','r')
trajs = []
if fr.mode=='r':
    filenames = fr.read().splitlines()
    for file in filenames:
        print(file)
        temp = np.load('%s'%file)
        trajs.append(temp)

lisa=[]
lisb=[]
lisc=[]
k=len(trajs)
print(k)

for i in range(0,84):
    if i < 28:
        lisa.append(i)
    if 28<=i<56:
        lisb.append(i)
    if 56<=i<84:
        lisc.append(i)

for i in range(84,192):
    if i < 120:
        lisa.append(i)
    if 120<=i<156:
        lisb.append(i)
    if 156<=i<192:
        lisc.append(i)

for i in range(192,228):
    if i < 204:
        lisa.append(i)
    if 204<=i<216:
        lisb.append(i)
    if 216 <= i <228:
        lisc.append(i)

#order O H C, 6+1 8+1 3
#252/3=84, C:3*4*3=12*3=36;H:9*4*3=36*3=108;O:7*4*3=28*3=84.
trajsa = []
trajsb = []
trajsc = []
for i in range(k):
    trajsa = trajs[i][:,lisa]
    trajsb = trajs[i][:,lisb]
    trajsc = trajs[i][:,lisc]
    np.save("outfile-short-a-C%02d"%(i),trajsa)
    np.save("outfile-short-b-C%02d"%(i),trajsb)
    np.save("outfile-short-c-C%02d"%(i),trajsc)

