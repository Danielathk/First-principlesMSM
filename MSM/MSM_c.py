#!/usr/bin/env pypy

'''
implied timescale analysis
'''

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
from msmbuilder.featurizer import AtomPairsFeaturizer
from msmbuilder.cluster import KCenters

from msmbuilder.io import load_trajs, save_trajs, save_generic
mp.use('Agg')

ntrajs=96
dim=2
nclusters=27
micro_trajs= []
for i in range(ntrajs):
    temp=np.load('micro_trajs-kmeans-%d-%d/%06d.npy'%(dim,nclusters,i))
    micro_trajs.append(temp)
micro_trajs = np.array(micro_trajs)

from msmbuilder.msm import MarkovStateModel, implied_timescales


lag_times=range(400,20000,400)
msm_timescales = implied_timescales(micro_trajs, lag_times, n_timescales=6,msm=MarkovStateModel(lag_time=250,reversible_type='transpose',ergodic_cutoff='off'))
np.savetxt('msm_timescales_1.txt',msm_timescales)

