#!/usr/bin/env pypy

'''
lumping of microstates to macrostate
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
lag_times=range(500,5000,500)
micro_trajs= []
for i in range(ntrajs):
    temp=np.load('micro_trajs-kmeans-%d-%d/%06d.npy'%(dim,nclusters,i))
    micro_trajs.append(temp)
micro_trajs = np.array(micro_trajs)

from msmbuilder.msm import MarkovStateModel, implied_timescales

print('building MSM')

msm = MarkovStateModel(lag_time=1000, reversible_type='transpose',ergodic_cutoff='off')
#msm = MarkovStateModel(lag_time=1000, reversible_type='mle')
msm.fit(micro_trajs)
print('finished MSM')
from msmbuilder.lumping import PCCAPlus
pcca = PCCAPlus.from_msm(msm, n_macrostates=12)
macro_trajs = pcca.transform(micro_trajs)
print('finished pcca+')

np.savetxt('pcca_2.txt', pcca.microstate_mapping_)
print('finished pcca+')
np.savetxt('membership_function_of_pccaplus.dat', pcca.chi_)
#np.savetxt('right_evs.adt', pcca.right_eigenvectors_)
np.savetxt('transition_matrix_20_1000.dat',msm.transmat_)



for i in range(ntrajs):
    np.savetxt('macro_trajs_%d'%i,macro_trajs[i],fmt='%2d')
