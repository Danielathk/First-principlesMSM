#!/usr/bin/env pypy

'''
perform GMRQ to select the number of k-centers
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

from sklearn.model_selection import KFold, RepeatedKFold

from sklearn.pipeline import Pipeline
from msmbuilder.cluster import NDGrid
from msmbuilder.msm import MarkovStateModel
import numpy as np

mp.use('Agg')

ntrajs=96
dim=2
nclusters=22
micro_trajs= []
for i in range(ntrajs):
    temp=np.load('micro_trajs-kmeans-%d-%d/%06d.npy'%(dim,nclusters,i))
    micro_trajs.append(temp)
micro_trajs = np.array(micro_trajs)


#n_states = [100, 300,500,700,1000,1500,2000, 2500, 3000]

lagtimes = [100, 200,300, 400,500, 600,700, 800,900, 1000]
cv = RepeatedKFold(n_splits=5, n_repeats=4, random_state=7)
results = []
model = Pipeline([
    ('msm', MarkovStateModel(n_timescales=20000, lag_time=10, reversible_type='transpose', verbose=False))
])

for lagtime in lagtimes:
    print("for #state = %d"%(lagtime))
    model = Pipeline([
    ('msm', MarkovStateModel(n_timescales=7, lag_time=lagtime, reversible_type='transpose', verbose=False))
]) 
    print(model)
    for fold, (train_index, test_index) in enumerate(cv.split(micro_trajs)):
        print('we are now dealing with fold%d'%(fold))
        train_data = [micro_trajs[i] for i in train_index]
        test_data = [micro_trajs[i] for i in test_index]

        # fit model with a subset of the data (training data).
        # then we'll score it on both this training data (which
        # will give an overly-rosy picture of its performance)
        # and on the test data.
        model.fit(train_data)
        train_score = model.score(train_data)
        test_score = model.score(test_data)
        
        print('training_score for this fold', train_score)
        print('testing score for this fold', test_score)

        results.append({
            'train_score': train_score,
            'test_score': test_score,
            'lagtimes': lagtime,
            'fold': fold})

import pandas as pd
results = pd.DataFrame(results)
results.head()

avgs = (results.groupby('lagtimes').aggregate(np.median).drop('fold', axis=1))

best_lagtime = avgs['test_score'].argmax()
best_score = avgs.loc[best_lagtime, 'test_score']
print(best_lagtime, "states gives the best score:", best_score)


#%matplotlib inline
from matplotlib import pyplot as plt

plt.scatter(results['lagtimes'], results['train_score'], c='b', lw=0, label=None)
plt.scatter(results['lagtimes'], results['test_score'],  c='r', lw=0, label=None)

plt.plot(avgs.index, avgs['test_score'], c='r', lw=2, label='Mean test')
plt.plot(avgs.index, avgs['train_score'], c='b', lw=2, label='Mean train')

plt.plot(best_lagtime, best_score, c='y', 
         marker='*', ms=20, label='{}'.format(best_lagtime))

plt.xscale('log')
plt.xlim((min(lagtimes)*.5, max(lagtimes)*5))
plt.ylabel('Generalized Matrix Rayleigh Quotient (Score)')
plt.xlabel('Lag Time')

plt.legend(loc='lower right', numpoints=1)
plt.tight_layout()
plt.savefig('GMRQ1.png', dpi=600)
plt.close()

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(12, 8.5))

# generate some random test data
#ll_data = [np.random.normal(0, std, 100) for std in range(6, 10)]
 #[100, 300,500,700,1000,1500,2000, 2500, 3000]
    
# plot violin plot
data_plot=[]
for j in [100, 200,300,400,500,600,700, 800, 900, 1000]:
    i=j//100-1
   # print(j)
    data_plot.append([])
    for k in range(20):
        data_plot[i].append(results['test_score'][20*i+k])
        
train_plot=[]
for j in [100, 200,300,400,500,600,700, 800, 900, 1000]:
    i=j//100-1
   # print(i)
    train_plot.append([])
    for k in range(20):
        train_plot[i].append(results['train_score'][20*i+k])

#axes[0].violinplot(data_plot,
#                   showmeans=False,
#                   showmedians=True)
#axes[0].violinplot(train_plot,
#                   showmeans=False,
 #                  showmedians=True)
#axes[0].set_title('violin plot')

# plot box plot

red_square = dict(markerfacecolor='r', marker='s')
trainplot=axes.boxplot(data_plot,patch_artist=True)
axes.boxplot(train_plot)
#axes.set_title('box plot')

# adding horizontal grid lines

axes.yaxis.grid(True)
axes.set_xticks([y+1 for y in range(len(data_plot))])
axes.set_xlabel('Lag time',size=22)
axes.set_ylabel('Generalized Matrix Rayleigh Quotient (Score)',size=22)
axes.set_ylim(5,6)


# add x-tick labels
plt.setp(axes, xticks=[y+1 for y in range(len(data_plot))],
         xticklabels=['100', '200','300', '400', '500', '600', '700', '800', '900','1000'])


#plt.show()
plt.savefig('GMRQ2.png', dpi=600)



