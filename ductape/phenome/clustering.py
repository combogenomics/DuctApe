#!/usr/bin/env python
"""
clustering

Phenome library

biolog data clustering functions
Many thanks to the scikits.learn team for the exhaustive documentation
"""
from itertools import product
from sklearn.cluster import KMeans, MeanShift, estimate_bandwidth
import numpy as np
import logging
import warnings

__author__ = "Marco Galardini"

logger = logging.getLogger('ductape.clustering')

def plotElbow(d, param_labels):
    from scipy.interpolate import interp1d
    import matplotlib.pyplot as plt
    
    logger.info('Saving elbow test plots')
    
    x_new = np.linspace(-1, max(d.keys()), 100)
    
    figidx = 1
    
    figsize = (len(d[d.keys()[0]][0])/2) + (len(d[d.keys()[0]][0])%2)  
    
    fig = plt.figure(figsize=(3.5*figsize, 8))
    fig.clf()
    for j in range(len(d[d.keys()[0]][0])):
        ax = fig.add_subplot(2, figsize, figidx)
        
        figidx += 1
        
        diffs = {}
        for i in d:
            diffs[i] = np.array([p[j] for p in d[i]]).mean()
        
        inter = interp1d(d.keys(), [diffs[i] for i in d], bounds_error=False,
                     kind='cubic')
        ax.plot(d.keys(),[diffs[i] for i in d],'o', x_new, inter(x_new),'-')
        ax.set_ylabel('Sum of squared errors')
        ax.set_xlabel('Num. clusters')
        ax.set_title(param_labels[j])
        
    fig.tight_layout()
    fig.savefig('elbow.png',dpi=300)

def plotClusters(X, labels, params=None, method='', prefix='clusters'):
    from ductape.common.utils import slice_it
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    
    X = np.array(X)
    labels = np.array(labels)
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)
    
    figidx = 1
    fig = plt.figure(1)
    fig.clf()
    for x, y in product(range(len(X[0])), repeat=2):
        ax = fig.add_subplot(len(X[0]),len(X[0]),figidx)
        
        if figidx%len(X[0]) == 1:
            if not params:
                ax.set_ylabel(x, rotation='horizontal')
            else:
                ax.set_ylabel(params[x], rotation='horizontal')
        if abs((len(X[0])*len(X[0])) % figidx - len(X[0])) <= len(X[0]):
            if not params:
                ax.set_xlabel(y)
            else:
                ax.set_xlabel(params[y])
        
        figidx += 1

        color = dict()
        j = 0
        for i in slice_it(range(255), cols=n_clusters_):
            color[j] = cm.RdYlGn(i[0])
            j += 1

        for k in range(n_clusters_):
            my_members = labels == k
            ax.plot(X[my_members, y], X[my_members, x], '.', color=color[k])

    fig.subplots_adjust(wspace=0, hspace=0)
    for ax in fig.axes:
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
    fig.suptitle('Clusters (%s, %s): %d' % (prefix, method, n_clusters_))
    fig.savefig('%s_%s.png'%(prefix,method),dpi=300)

def mean(X, save_fig=False, params_labels=None, prefix='clusters'):
    '''
    Compute clustering with MeanShift
    '''
    logger.debug('Calculating MeanShift clusters using %d parameters'%len(X[0]))
    
    X = np.array( X )
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bandwidth = estimate_bandwidth(X, quantile=0.2)
    
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(X)
        
    labels = ms.labels_
    
    if save_fig:
        plotClusters(X, ms, method='mean', prefix=prefix,
                     params=params_labels)
    
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)
    
    logger.debug('Found %d clusters with MeanShift algorithm'%n_clusters_)
    
    return labels

def _kmeans(X, n_clusters):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        
        try:
            k_means = KMeans(init='random', k=n_clusters, n_init=100,
                             max_iter=1000)
        except:
            k_means = KMeans(init='random', n_clusters=n_clusters, n_init=100,
                             max_iter=1000)
        k_means.fit(X)
        
    return k_means

def getSseKmeans(k_means, X):
    '''
    Returns a list of sums of squared errors
    '''
    import math
    
    sse=[]
    for j in range(len(k_means.labels_)):
        c = k_means.labels_[j]
        ref = k_means.cluster_centers_[c]
        dist = []
        for k in range(len(ref)):
            dist.append(math.pow(X[j][k] - ref[k], 2))
        sse.append(dist)
        
    return sse

def kmeans(X, n_clusters=10, save_fig=False, params_labels=None, prefix='clusters'):
    '''
    Compute clustering with KMeans
    '''
    logger.debug('Calculating KMean clusters using %d parameters'%len(X[0]))
    
    X = np.array( X )
    
    k_means = _kmeans(X, n_clusters)
    
    labels = k_means.labels_
    
    if save_fig:
        plotClusters(X, k_means, method='kmeans', prefix=prefix,
                     params=params_labels)
    
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)
    
    logger.debug('Found %d clusters with KMeans algorithm'%n_clusters_)
    
    return labels
    