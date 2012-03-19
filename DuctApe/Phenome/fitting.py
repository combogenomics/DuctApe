#!/usr/bin/env python
"""
fitting

Phenome library

biolog data fitting functions
"""
import numpy as np

def gompertz(x, A, u, d, y0):
    '''
    Gompertz growth model
    Taken from: "Modeling of the bacterial growth curve."
                (Zwietering et al., 1990)
                PMID: 16348228
    '''
    y = (A * np.exp( -np.exp( (((u * np.e)/A) * (d - x)) + 1 ) ) ) + y0
    return y

def logistic(x, A, u, d, y0):
    '''
    Logistic growth model
    Taken from: "Modeling of the bacterial growth curve."
                (Zwietering et al., 1990)
                PMID: 16348228
    '''
    y = (A / (1 + np.exp( ( ((4 * u)/A) * (d - x) ) + 2 ))) + y0
    return y

def getFlex(x, y):
    '''
    Given two axes (with the same length!) returns a guess of the flex point
    '''
    if len(x) != len(y):
        raise ValueError('Axes have different sizes (x: %d, y: %d)'%(len(x),len(y)))
    
    diffs = []
    indexes = range(len(x))
    
    for i in indexes:
        if i+1 not in indexes:
            continue
        diffs.append(y[i+1] - y[i])
    diffs = np.array( diffs )
    
    flex = x[-1]
    for i in indexes:
        if i+1 not in indexes:
            continue
        if (y[i+1] - y[i]) > (diffs.mean() + (diffs.std())):
            flex = x[i]
            break
    
    return flex