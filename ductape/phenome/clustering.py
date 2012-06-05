#!/usr/bin/env python
"""
clustering

Phenome library

biolog data clustering functions
Many thanks to the scikits.learn team for the exhaustive documentation
"""
from itertools import cycle, product
from sklearn.cluster import KMeans, MeanShift, estimate_bandwidth
import numpy as np
import logging
import warnings
import matplotlib.pyplot as plt

__author__ = "Marco Galardini"

logger = logging.getLogger('ductape.clustering')

def plotClusters(X, clust, params=None, method='', prefix='clusters'):
    labels = clust.labels_
    cluster_centers = clust.cluster_centers_

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

        colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
        for k, col in zip(range(n_clusters_), colors):
            if k > len(cluster_centers) - 1:continue
            my_members = labels == k
            cluster_center = cluster_centers[k]
            ax.plot(X[my_members, y], X[my_members, x], col + '.')
            ax.plot(cluster_center[y], cluster_center[x], 'o', markerfacecolor=col,
                                            markeredgecolor='k', markersize=2)

    fig.subplots_adjust(wspace=0, hspace=0)
    for ax in fig.axes:
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
    fig.suptitle('Clusters: %d' % n_clusters_)
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
        plotClusters(X, ms, method='mean', prefix=prefix)
    
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)
    
    logger.debug('Found %d clusters with MeanShift algorithm'%n_clusters_)
    
    return labels

def kmeans(X, n_clusters=10, save_fig=False, params_labels=None, prefix='clusters'):
    '''
    Compute clustering with KMeans
    '''
    logger.debug('Calculating KMean clusters using %d parameters'%len(X[0]))
    
    X = np.array( X )
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        
        k_means = KMeans(init='random', k=n_clusters, n_init=100, max_iter=1000)
        k_means.fit(X)
    
    labels = k_means.labels_
    
    if save_fig:
        plotClusters(X, k_means, method='kmeans', prefix=prefix)
    
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)
    
    logger.debug('Found %d clusters with KMeans algorithm'%n_clusters_)
    
    return labels
    