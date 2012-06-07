#!/usr/bin/env python
"""
fitting

Phenome library

biolog data fitting functions
"""
from ductape.common.utils import compress, smooth
from scipy.optimize.minpack import curve_fit
import numpy as np
import logging
# No country for warnings
import scipy as sp
sp.seterr(all='ignore')
#

__author__ = "Marco Galardini"

logger = logging.getLogger('ductape.fitting')

def logistic(x, A, u, d, v, y0):
    '''
    Logistic growth model
    Taken from: "Modeling of the bacterial growth curve."
                (Zwietering et al., 1990)
                PMID: 16348228
    '''
    y = (A / (1 + np.exp( ( ((4 * u)/A) * (d - x) ) + 2 ))) + y0
    return y

def gompertz(x, A, u, d, v, y0):
    '''
    Gompertz growth model
    Taken from: "Modeling of the bacterial growth curve."
                (Zwietering et al., 1990)
                PMID: 16348228
    '''
    y = (A * np.exp( -np.exp( (((u * np.e)/A) * (d - x)) + 1 ) ) ) + y0
    return y

def richards(x, A, u, d, v, y0):
    '''
    Richards growth model
    (equivalent to Stannard)
    Taken from: "Modeling of the bacterial growth curve."
                (Zwietering et al., 1990)
                PMID: 16348228
    '''
    y = (A * pow(1 + (v + (np.exp(1 + v) * np.exp( (u/A) * (1 + v) * (1 + (1/v)) * (d - x) ) ) ),-(1/v))) + y0
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

def getPlateau(x, y):
    '''
    Given two axes (with the same length!) returns a guess of the plateau point
    '''
    if len(x) != len(y):
        raise ValueError('Axes have different sizes (x: %d, y: %d)'%(len(x),len(y)))
    
    ymax = y.max()
    
    diffs = []
    indexes = range(len(y))
    
    for i in indexes:
        if i+1 not in indexes:
            continue
        diffs.append(y[i+1] - y[i])
    diffs = np.array( diffs )
    
    ymax = y[-1]
    for i in indexes:
        if y[i] > (ymax - diffs.std()) and y[i] < (ymax + diffs.std()):
            ymax = y[i]
            break
    
    return ymax
    

def rect(x, a, y0):
    '''
    yep, that's a rect!
    '''
    y = (a * x) + y0
    return y

def fitData(xdata, ydata):
    '''
    Fits the provided data to the first working function
    (first Gompertz, then Logistic, then Richards)

    Returns a tuple with plateau, slope, lag, y0 and model used
    If no fitting was possible all values are None

    Please note that the plateau may be reached outside the final time point
    '''
    retries = 3
    while retries > 0:
        params = [None, None, None, None, None]
        model = ''
        # Initial guesses for the output parameters
        p0 = [getPlateau(xdata, ydata), 4.0, getFlex(xdata, ydata), 0.1, 0]
        if retries == 1:
            p0[2] = 0
        try:
            params, pcov = curve_fit(gompertz, xdata, ydata, p0 = p0)
            model = 'gompertz'
            break
        except:
            #logger.debug('Gompertz fit failed')
            try:
                params, pcov = curve_fit(logistic, xdata, ydata, p0 = p0)
                model = 'logistic'
                break
            except:
                #logger.debug('Logistic fit failed')
                try:
                    params, pcov = curve_fit(richards, xdata, ydata, p0 = p0)
                    model = 'richards'
                    break
                except:
                    #logger.debug('Richards fit failed')
                    retries -= 1
                    #logger.debug('%d retries left'%retries)
                    # Compress again the data
                    ydata = np.array(compress(ydata, span=4))
                    ydata = np.array(smooth(ydata, window_len = len(ydata)/2, 
                              window = 'blackman'))
                    xdata = np.array(compress(xdata, span=4))
                    #
                    params = [None, None, None, None, None]
    
    return params, model