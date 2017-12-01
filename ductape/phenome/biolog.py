#!/usr/bin/env python
"""
Biolog

Phenome library

Classes to handle Biolog data
"""
from ductape import __email__
from ductape.common.commonthread import CommonThread
from ductape.common.utils import smooth, compress
from matplotlib import cm
from matplotlib import colors
import Queue
import csv
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
# No country for warnings
np.seterr(all='ignore')
#

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('ductape.biolog')

################################################################################
# Classes

class Well(object):
    '''
    Class Well
    Contains signals for a particular plate/well
    '''
    def __init__(self, plate_id, well_id):
        self.plate_id = plate_id
        self.well_id = well_id.replace(' ','')
        self.signals = {}
        self.smoothed = False
        self.compressed = False
        
        # Parameters
        self.max = None
        self.min = None
        self.height = None
        self.plateau = None
        self.slope = None
        self.lag = None
        self.area = None
        self.v = None
        self.y0 = None
        
        # Parameters list
        self.params = ['max', 'min', 'height',
                       'plateau', 'slope', 'lag',
                       'area', 'v', 'y0']
        
        # Fitting model used
        self.model = None
        
        # Source of the parameters (i.e. the program)
        self.source = None
        
        # Relative activity index
        self.activity = None
        
        self.otherparams = ['activity', 'model', 'source']
        
        # Additional info added by well parents
        self.replica = None
        self.strain = None
        self.zero = False

    def getHeader(self):
        '''
        Get the header of the returned value from str(Well)
        '''
        return '\t'.join( ['Plate', 'Well', 'Strain',
                           'Replica', 'Activity', 'Min',
                           'Max', 'Height', 'Plateau',
                           'Slope', 'Lag', 'Area', 'Source'] )
        
    def __str__(self):
        '''
        Note: the header can be retrieved by calling getHeader
        '''
        return '\t'.join( [self.plate_id, self.well_id, self.strain] +
                          [str(x) for x in [self.replica,
                                            self.activity,
                                            self.min,
                                            self.max,
                                            self.height,
                                            self.plateau,
                                            self.slope,
                                            self.lag,
                                            self.area,
                                            self.source]] )

    def addSignal(self,time,signal):
        self.signals[time] = signal
        
    def fillMissing(self, times):
        '''
        Given a times list, fills the missing values with 0
        '''
        for hour in sorted(times):
            if hour not in self.signals:
                if sorted(times).index(hour) - 1 >= 0:
                    self.signals[hour] = self.signals[
                                          sorted(times)[
                                                sorted(times).index(hour) - 1]]
                else:
                    self.signals[hour] = np.NAN
    
    def getMax(self):
        '''
        Maximum signal
        '''
        return max(self.signals.values())
    
    def getMin(self):
        '''
        Minimum signal
        '''
        return min(self.signals.values())
    
    def getMaxTime(self):
        '''
        Maximum signal
        '''
        return max(self.signals.keys())
    
    def getMinTime(self):
        '''
        Minimum signal
        '''
        return min(self.signals.keys())
    
    def smooth(self, window_len = 11, window_type = 'hanning',
               forceZero = True):
        '''
        Apply a smoothing algorithm to the signals
        Really useful for clearer plots and other features
        Available windows: 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
        '''
        if not self.smoothed:
            signals = [self.signals[hour] for hour in sorted(self.signals.keys())]
            
            # If there are not enough signals, do not smooth
            if len(signals) <= 3*11:
                logger.debug('Too few time points for %s %s (%d): no smoothing'%
                            (self.plate_id, self.well_id, len(signals)))
                return

            smoothed = smooth(signals, window_len = window_len, 
                              window = window_type)
            
            for idx in range(len(self.signals)):
                hour = sorted(self.signals.keys())[idx]
                if smoothed[idx] < 0 and forceZero:
                    self.signals[hour] = 0.1
                else:
                    self.signals[hour] = smoothed[idx]
            
            self.smoothed = True
        else:
            logger.warning('Plate %s, Well %s was already smoothed'%
                          (self.plate_id, self.well_id))
            
    def compress(self, span = 3):
        '''
        Reduce the amount of signals
        This function should be called BEFORE smooth
        '''
        if self.smoothed:
            logger.warning('Plate %s, Well %s should be smoothed AFTER compression'%
                          (self.plate_id, self.well_id))
        
        if not self.compressed:
            times = compress( sorted(self.signals.keys()), span = span)
            
            # If there are not enough time points, do not compress
            if len(times) <= 3*11:
                logger.debug('Too few time points for %s %s (%d): no compress'%
                            (self.plate_id, self.well_id, len(times)))
                return
            
            toremove = [t for t in self.signals.keys() if t not in times]
            for t in toremove:
                del self.signals[t]
            
            self.compressed = True
        else:
            logger.warning('Plate %s, Well %s was already compressed'%
                          (self.plate_id, self.well_id))
    
    def calculateParams(self,
                            noCompress = False, noSmooth = False):
        '''
        Populates the parameters values for the experiment
        By default compression and smoothing are applied to save some time
        '''
        from scipy.integrate import trapz
        from ductape.phenome.fitting import fitData, getFlex, getPlateau
       
        if not self.compressed:
            self.compress()
        if not self.smoothed:
            self.smooth(window_len=11, window_type='blackman')
            
        # Let's start with the easy ones!
        self.max = self.getMax()
        
        self.min = self.getMin()
        
        self.height = np.array( self.signals.values() ).mean()
        
        # Let's go with the function fitting
        xdata = np.array( [x for x in sorted(self.signals.keys())] )
        ydata = np.array( [self.signals[x] for x in xdata] )
        (self.plateau, self.slope, self.lag, v, y0), self.model = fitData(xdata, ydata)
        
        # May be needed for debugging purposes
        # or to plot some fitting data
        self.v = v
        self.y0 = y0
        
        # Trapezoid integration for area calculation
        self.area = trapz(y = ydata, x = xdata)
        
        self.source = "DuctApe"
        
        # If any of the values are null generate them by hand
        if not y0:
            self.plateau = 0
            self.slope = 0
            self.lag = 0
            self.v = 0
            self.y0 = 0
            return
        
        # Check the fitting parameters
        if self.slope < 0 or self.slope > ydata.max() or self.plateau < 0:
            self.plateau = 0
            self.lag = 0
            self.slope = 0.0
        else:
            if self.lag >= 0:
                y0 = - (self.lag * self.slope)
            else:
                y0 = 0
                self.lag = 0
            xplateau = (self.plateau - y0) / self.slope
            if xplateau > xdata.max():
                self.plateau = getPlateau(xdata, ydata)
                self.lag = getFlex(xdata, ydata)
                
                xplat = list(xdata)[list(ydata).index(self.plateau)]
                ylag = list(ydata)[len(xdata) - list(xdata)[::-1].index(self.lag) - 1]
                
                self.slope = np.sqrt(pow((xplat - self.lag), 2) +
                                     pow((self.plateau - ylag), 2))
                
                if self.slope < 0 or self.slope > ydata.max() or self.plateau < 0:
                    self.plateau = 0
                    self.lag = 0
                    self.slope = 0.0
                    
    def purgeNan(self):
        '''
        Scan the parameters list, force Nan values to zero
        Works also for missing parameters (None)
        '''
        import math
        
        for param in self.params:
            if getattr(self, param) is None:
                logger.debug('Well %s %s %s, parameter %s had a None value,'
                           %(self.plate_id, self.well_id, self.strain, param)+
                           ' forced to zero')
                setattr(self, param, 0)
            if math.isnan( float( getattr(self, param) ) ):
                logger.debug('Well %s %s %s, parameter %s had a NaN value,'
                           %(self.plate_id, self.well_id, self.strain, param)+
                           ' forced to zero')
                setattr(self, param, 0)
                
    def isParams(self):
        '''
        Do we have at least one parameter calculated?
        '''
        params = filter(lambda x: x!=None,
                      [self.max,
                       self.min,
                       self.height,
                       self.plateau,
                       self.slope,
                       self.lag,
                       self.area,
                       self.v,
                       self.y0])
        
        if len(params) == 0:
            return False
        else:
            return True
        
    def hasMissingParams(self):
        return None in set([self.max,
                           self.min,
                           self.height,
                           self.plateau,
                           self.slope,
                           self.lag,
                           self.area,
                           self.v,
                           self.y0])

class SinglePlate(object):
    '''
    Class SinglePlate
    Contains informations about a particular plate
    '''
    def __init__(self):
        self.plate_id = None
        self.strainType = None
        self.sample = None
        self.strainName = None
        self.strainNumber = None
        self.other = None
        
        # Added by the system to avoid confusion
        self.strain = None
        
        # Replica management
        self.replica = None
        
        # Raw data --> well_id -> Well objects
        self.data = {}
        # Used internally, index in Hours row --> well_id
        self._idx = {}
        
        # This plate was zero-subtracted?
        self.zero = False
        
    def getMax(self):
        '''
        Maximum signal for the entire plate
        '''
        return max( [self.data[well].getMax() for well in self.data] )
    
    def calculateParams(self):
        '''
        Iterate over each well: calculate parameters and return the
        A generator is returned
        '''
        for well_id, well in self.data.iteritems():
            if not well.isParams():
                well.calculateParams()
            yield True
    
    def getWells(self):
        '''
        Generator to get the single wells
        '''
        for well_id, well in self.data.iteritems():
            well.replica = self.replica
            well.strain = self.strain
            yield well
            
    def addWell(self,well):
        '''
        Add a Well object
        If it is already present, a warning will pop-up but it will be replaced
        '''
        if well.well_id in self.data:
            logger.warning('Replacing an already existing well (%s, %s, %s, %s)'
                            %(self.plate_id, well.well_id, self.strain,
                              self.replica))
        
        self.data[well.well_id] = well
            
class Plate(object):
    '''
    Class Plate
    Contains all the distinct strains data for a particular plate
    There can be more than one replica for each strain
    '''
    def __init__(self, plate_id, plate_name = '', smooth = True, window = 11,
                 alpha = 0.5, compress = 0,
                 linewidth = 2, maxsig = None):
        self.plate_id = plate_id
        self.plate_name = ' '.join( (plate_id, plate_name) ).rstrip()
        self.strains = {}
        self.colors = {}
        self.wellNames = {}
        
        self.times = None
        self.wells = None
        
        self.smooth = bool(smooth)
        self.window = int(window)
        self.compress = int(compress)
        
        self.alpha = float(alpha)
        self.linewidth = float(linewidth)
        
        self.maxsignal = maxsig
        
        # Figure(s)
        self.figure = None
        self.heatfig = None
        self.legend = None
        
        self._figidx = 1
    
    def getMax(self):
        return max([plate.getMax() 
                    for strain, plates in self.strains.iteritems()
                    for plate in plates])
                    
    def getMaxActivity(self):
        return max([w.activity for w in self.getWells()])
        
    def calculateParams(self):
        '''
        Iterate over the single wells and calculate the parameters
        '''
        for strain, plates in self.strains.iteritems():
            for plate in plates:
                plate.calculateParams()
            yield True
                    
    def getWells(self):
        '''
        Generator to get the single wells
        '''
        for strain, plates in self.strains.iteritems():
            for plate in plates:
                for well in plate.getWells():
                    yield well
    
    def addWellTitles(self, dWell):
        '''
        Add all the titles of the current plate in the form of a dictionary
        d[well_id] = title
        '''
        self.wellNames = dWell
    
    def _bracketing(self, signals):
        '''
        Given a set of signals returns two sets with maximum and minimum value
        '''
        if len(signals) == 2:
            return signals
        
        maxout = []
        minout = []
        
        for i in range(len(signals[0])):
            sigs = [x[i] for x in signals]
            maxout.append(max(sigs))
            minout.append(min(sigs))
            
        return maxout, minout
    
    def _plot(self, well_id, dWell, ax):
        '''
        Smooths the signal and plots it
        If there are more than one exp for a strain, the intersection is plotted
        '''
        for strain,signals in dWell.iteritems():
            if len(signals) > 1:
                # Intersect!
                maxsig,minsig = self._bracketing(signals)
                if self.smooth:
                    maxsig = smooth(maxsig, window_len=self.window)
                    minsig = smooth(minsig, window_len=self.window)
                ax.fill_between(self.times, maxsig, minsig,
                                 color=self.colors[strain],
                                 linewidth=self.linewidth,
                                 alpha=self.alpha,
                                 rasterized=True)
            elif len(signals) == 0:
                continue
            else:
                # Single plot!
                if self.smooth:
                    signal = smooth(signals[0], window_len=self.window)
                else:
                    signal = signals[0]
                ax.plot(self.times, signal, color=self.colors[strain],
                        linewidth=self.linewidth, rasterized=True)
                    
        ax.set_ylim(0,self.maxsignal)
        
    def fixFigure(self):
        '''
        Fix some parameters of the Big picture
        '''
        self.figure.subplots_adjust(wspace=0, hspace=0)
        for ax in self.figure.axes:
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
        self.figure.suptitle(self.plate_name)
    
    def preparePlot(self):
        '''
        Prepare a series of plots
        '''
        # Check colors
        for strain in self.strains:
            if strain not in self.colors:
                logger.error('Color code for strain %s is missing!'%strain)
                return
                
        # Check time concordance
        times = []
        for strain, plates in self.strains.iteritems():
            for plate in plates:
                for well_id, data in plate.data.iteritems():
                    for time in data.signals.keys():
                        if time not in times:
                            times.append(time)
                    break
                
        times.sort()
        if self.compress != 0:
            times = compress(times, self.compress)
        self.times = times
        
        # Get also each plate/well pair
        wells = []
        for strain, plates in self.strains.iteritems():
            for plate in plates:
                for well_id, data in plate.data.iteritems():
                    data.fillMissing(self.times)
                    if well_id not in wells:
                        wells.append(well_id)
        wells.sort()
        self.wells = wells
        
    def _prepareSignals(self, well_id):
        '''
        Prepares the signals for a specific well
        '''
        strain_signals = {}
        
        for strain, plates in self.strains.iteritems():
            strain_signals[strain] = []
            for plate in plates:
                try:
                    strain_signals[strain].append( 
                            [plate.data[well_id].signals[hour]
                             for hour in self.times])
                except:
                    logger.debug('Something missing: %s, %s, %f'%(
                                    strain, well_id, hour))
            
        return strain_signals
    
    def plotAll(self):
        # Preparatory steps
        if not self.times and not self.wells:
            self.preparePlot()
        
        self._figidx = 1
        
        # Cycle over each well
        for well_id in self.wells:
            strain_signals = self._prepareSignals(well_id)
            
            if not self.figure:
                self.figure = plt.figure()
            ax = self.figure.add_subplot(8,12,self._figidx)
            
            # Smooth & Plot
            self._plot(well_id, strain_signals, ax)
            
            if self._figidx%12 == 1:
                ax.set_ylabel(well_id[0], rotation='horizontal')
            if (96 - self._figidx < 12) and abs(96 % self._figidx - 12) <= 12:
                ax.set_xlabel(str(abs(96 % self._figidx - 12)))
            self._figidx += 1
            
            yield well_id, self.wells.index(well_id)
        
        self.fixFigure()
    
    def arrayReshape(self, acts, strains):
        # Array reshape (tricky!)
        # 1- Finger crossed for a perfect square
        square = np.sqrt( len(acts) )
        if square.is_integer():
            # Yep
            acts = acts.reshape(int(square), int(square))
        else:
            # Oh snap!
            cols = int(square)
            rows = len(acts) / float(cols)
            if rows.is_integer():
                acts = acts.reshape(int(rows), cols)
            else:
                # Some fake signals will be added
                rows = int(rows) + 1
                new = (cols * rows) - len(acts)
                acts = np.array( acts.tolist() + [np.nan for i in range(new)] )
                acts = acts.reshape(rows, cols)
        
        return acts
    
    def plotActivity(self, strains=[], maxAct=9):
        '''
        Generator:
        Plots the activity for each strain as heatmaps
        A strains subset can be provided, otherwise all strains are plotted 
        in alphabetival order 
        '''
        # Reality check on provided strains
        if len(strains) > 0:
            strains = set(strains)
            unknown = set(strains).difference(self.strains.keys())
            if len( unknown ) > 0:
                logger.warning('Unknown strain(s) were provided (%s)'%
                                 ' '.join(unknown))
                # Fall back
                strains = sorted(self.strains.keys())
        else:
            strains = sorted(self.strains.keys())
        
        # Preparatory steps
        if not self.times and not self.wells:
            self.preparePlot()
        
        self._figidx = 1
        
        # Cycle over each well
        for well_id in self.wells:            
            if not self.heatfig:
                self.heatfig = plt.figure()
            ax = self.heatfig.add_subplot(8,12,self._figidx)
            
            # Plot activity matrix (AKA the almighty heatmap)
            # Array preparation
            acts = np.array([self.strains[strain][0].data[well_id].activity
                         for strain in strains])
            acts = self.arrayReshape(acts, strains)
                    
            ax.matshow(acts, cmap=cm.RdYlGn, vmin=0, vmax=maxAct)
            
            if self._figidx%12 == 1:
                ax.set_ylabel(well_id[0], rotation='horizontal')
            if (96 - self._figidx < 12) and abs(96 % self._figidx - 12) <= 12:
                ax.set_xlabel(str(abs(96 % self._figidx - 12)))
            self._figidx += 1
            
            yield well_id, self.wells.index(well_id)
        
        self.heatfig.subplots_adjust(wspace=0, hspace=0)
        for ax in self.heatfig.axes:
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
        self.heatfig.suptitle(self.plate_name)
    
    def plotLegend(self, plate_id, strains=[]):
        '''
        Generate a plot with the position of the strains in the activity plots
        and the color reference
        '''
        # Reality check on provided strains
        if len(strains) > 0:
            strains = set(strains)
            unknown = set(strains).difference(self.strains.keys())
            if len( unknown ) > 0:
                logger.warning('Unknown strain(s) were provided (%s)'%
                                 ' '.join(unknown))
                # Fall back
                strains = sorted(self.strains.keys())
        else:
            strains = sorted(self.strains.keys())
        
        array = np.array([strain for strain in strains])
        array = self.arrayReshape(array, strains)
        
        self.legend = plt.figure()
        ax = self.legend.add_subplot(111)
        
        # Show the array
        acts = np.array([0 for strain in strains])
        acts = self.arrayReshape(acts, strains)
        
        # Setup the legend
        ax.matshow(acts, cmap='Greys', vmin=0, vmax=1)
        for i in range(acts.shape[0] - 1):
            ax.axhline(y=i+0.5, color='black')
        
        for i in range(acts.shape[0] - 1):
            ax.axvline(x=i+0.5, color='black')
        
        # Plot the strains names in the proper order
        for i in range(len(array)):
            for j in range(len(array[i])):
                if array[i,j] not in self.colors:
                    continue
                ax.text(j, i, array[i,j], color=self.colors[array[i,j]],
                        fontsize='x-large', fontweight='bold', ha='center')
        
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Strains color codes and order in heatmap (%s)'%plate_id)
    
    def plotWell(self, well_id, fig=None):
        '''
        Generates and returns a single well as a figure
        '''
        # Preparatory steps
        if not self.times and not self.wells:
            self.preparePlot()
        
        strain_signals = self._prepareSignals(well_id)
        
        # Temporary increase in the line width
        self.linewidth += 3 
        
        # Figure creation
        if not fig:
            fig = plt.figure()
        ax = fig.add_subplot(111)
        
        self._plot(well_id, strain_signals, ax)
        
        if well_id in self.wellNames:
            ax.set_title(' '.join( [well_id, self.wellNames[well_id]]))
        ax.set_xlabel('Hour')
        ax.set_ylabel('Signal')
        
        self.linewidth -= 3
        
        return fig
        
    def setColor(self, strain, color):
        '''
        Set the color code of a strain
        '''
        self.colors[strain] = color
    
    def addData(self, strain, data):
        '''
        Add a SinglePlate object regarding a particular strain
        A check on the plate_id is performed!
        '''
        if data.plate_id != self.plate_id:
            logger.error('Expecting %s, got %s'%(self.plate_id,data.plate_id))
            return False
        
        if strain not in self.strains:
            self.strains[strain] = []

        self.strains[strain].append(data)
        
        # This should not be needed anymore...
        # Replica should be handled from the DB
        #data.replica = len(self.strains[strain])
        # Keeping the scarf tissue code just in case something breaks
        
        return True

class Experiment(object):
    '''
    Class Experiment
    Contains all the data (including replicas) for a distinct biolog experiment
    Can perform clusterization, replica management
    
    Input: plates [Plate object]
    '''
    class well(object):
        pass
    
    def __init__(self, exp_id='', name='', plates=[], zero=False,
                 category = {}, categorder = [], zeroPlates=set()):
        self.exp_id = exp_id
        self.name = name
        
        self.zero = zero
        self.category = category
        self.categorder = categorder
        self.zeroPlates = zeroPlates
        
        self.plates = {}
        for plate in plates:
            if not self._addPlate(plate):
                self.plates = {}
                break
        
        self.maxParams = {}
        
        self.experiment = {}
        self.sumexp = {}
        self._organize()
        
        # Allowed policies for purging of replicas
        self.policies = ['keep-min', 'keep-max',
                         'keep-min-one', 'keep-max-one',
                         'replica']
        
        self.purged = False
        
        self.discarded = set()
    
    def _addPlate(self, plate):
        if plate.plate_id not in self.plates:
            self.plates[plate.plate_id] = plate
            return True
        else:
            logger.warning('Plate %s already present! Please retry...'%
                                    plate.plate_id)
            return False
        
    def _organize(self):
        '''
        Organize the whole experiment in a dictionary-based structure
        '''
        for w in self.getWells(params=False):
            if w.plate_id not in self.experiment:
                self.experiment[w.plate_id] = {}
                self.sumexp[w.plate_id] = {}
            if w.well_id not in self.experiment[w.plate_id]:
                self.experiment[w.plate_id][w.well_id] = {}
                self.sumexp[w.plate_id][w.well_id] = {}
            if w.strain not in self.experiment[w.plate_id][w.well_id]:
                self.experiment[w.plate_id][w.well_id][w.strain] = {}
                
                fakeWell = Well(w.plate_id, w.well_id)
                fakeWell.strain = w.strain 
                self.sumexp[w.plate_id][w.well_id][w.strain] = fakeWell
            
            self.experiment[w.plate_id][w.well_id][w.strain][w.replica] = w
            
        # For each well, if each replica has an activity, add its value
        # to the averaged dictionary
        for pid in self.sumexp:
            for wid in self.sumexp[pid]:
                for org in self.sumexp[pid][wid]:
                    reps = self.experiment[pid][wid][org].keys()
                    
                    # Keep track of all mean parameters
                    for param in Well('phony', 'phony').params + ['activity']:
                        act = []
                        for r in reps:
                            if getattr(self.experiment[pid][wid][org][r], param) is None:
                                break
                            act.append(getattr(self.experiment[pid][wid][org][r], param))
                        
                        if len(act) > 0:
                            mean = np.array(act).mean()
                            setattr(self.sumexp[pid][wid][org], param, mean)
    
    def getMax(self):
        '''
        Get the maximum signal value of the whole experiment
        '''
        return max([plate.getMax()
                    for plate_id, Plate in self.plates.iteritems()
                    for strain, plates in Plate.strains.iteritems()
                    for plate in plates])
        
    def calculateParams(self):
        '''
        Generator to the single well parameters for clustering
        '''
        for plate_id in self.plates:
            Plate = self.plates[plate_id]
            for res in Plate.calculateParams():
                yield True
    
    def getRandomZeroWells(self, howmany=10, activity=None):
        '''
        Generator to random wells (zero subtracted)
        If activity is set, only those wells w/ desired AV will be retrieved
        '''
        import random
        
        # Check if activity makes sense
        if activity is not None:
            if activity not in self.getDistinctActivity():
                raise ValueError('Provided AV value not found in the wells pool')
        
        wells = []
        
        while howmany != 0:
            while True:
                while True:
                    plate_id = random.choice( list(self.experiment.keys()) )
                    if plate_id in self.zeroPlates:
                        break
                well_id = random.choice( list(self.experiment[plate_id].keys()) )
                strain = random.choice( 
                            list(self.experiment[plate_id][well_id].keys()) )
                replica = random.choice( 
                            list(self.experiment[plate_id][well_id][strain].keys()))
                
                w = self.experiment[plate_id][well_id][strain][replica]
                if activity is not None:
                    if w.activity == activity:
                        wells.append(w)
                    else:continue
                else:
                    wells.append(w)
                howmany -= 1
                break
        
        return wells
    
    def getRandomNoZeroWells(self, howmany=10, activity=None):
        '''
        Generator to random wells (no zero subtracted)
        If activity is set, only those wells w/ desired AV will be retrieved
        '''
        import random
        
        # Check if activity makes sense
        if activity is not None:
            if activity not in self.getDistinctActivity():
                raise ValueError('Provided AV value not found in the wells pool')
        
        wells = []
        
        while howmany != 0:
            while True:
                while True:
                    plate_id = random.choice( list(self.experiment.keys()) )
                    if plate_id not in self.zeroPlates:
                        break
                well_id = random.choice( list(self.experiment[plate_id].keys()) )
                strain = random.choice( 
                            list(self.experiment[plate_id][well_id].keys()) )
                replica = random.choice( 
                            list(self.experiment[plate_id][well_id][strain].keys()))
                
                w = self.experiment[plate_id][well_id][strain][replica]
                if activity is not None:
                    if w.activity == activity:
                        wells.append(w)
                    else:continue
                else:
                    wells.append(w)
                howmany -= 1
                break
        
        return wells

    def getZeroWells(self, params=True):
        '''
        Generator to get the zero-subtracted single wells
        if params is set to False, it just gives you the wells,
        otherwise it calculates them
        '''
        for well in self.getWells(params):
            if well.plate_id in self.zeroPlates:
                yield well
                
    def getNoZeroWells(self, params=True):
        '''
        Generator to get the nonzero-subtracted single wells
        if params is set to False, it just gives you the wells,
        otherwise it calculates them
        '''
        for well in self.getWells(params):
            if well.plate_id not in self.zeroPlates:
                yield well
                
    def getCategoryWells(self, params=True):
        '''
        Generator to get ('category', [wells])
        if params is set to False, it just gives you the wells,
        otherwise it calculates them
        '''
        for categ in self.categorder:
            plates = self.category[categ]
            wells = []
            for well in self.getWells(params):
                if well.plate_id in plates:
                    wells.append(well)
            yield (categ, wells)
    
    def getWells(self, params=True):
        '''
        Generator to get the single wells
        if params is set to False, it just gives you the wells,
        otherwise it calculates them
        '''
        for plate_id in self.plates:
            Plate = self.plates[plate_id]
            for well in Plate.getWells():
                if not well.isParams() and params:
                    well.calculateParams()
                
                yield well
    
    def setNoActivity(self):
        '''
        All the wells have no activity!
        '''
        # TODO: distinguish between zero and nonzero
        for plate_id in self.plates:
            Plate = self.plates[plate_id]
            for strain, plates in Plate.strains.iteritems():
                for plate in plates:
                    for wid, well in plate.data.iteritems():
                        well.activity = 0
    
    def getPurgedWells(self):
        '''
        Generator to get the single purged wells
        '''
        for plate in self.experiment:
            for well in self.experiment[plate]:
                for strain in self.experiment[plate][well]:  
                    reps = self.experiment[plate][well][strain].values()
                    for w in reps:
                        yield w
                        
    def getAverageWells(self, org_id=None):
        '''
        Generator to the single average wells
        if org_id is provided, only the wells from that ID is provided
        (plate IDs, well IDs and organism IDs are sorted)
        '''
        for plate in sorted(self.sumexp.keys()):
            for well in sorted(self.sumexp[plate].keys()):
                if org_id:
                    if org_id not in self.sumexp[plate][well]:continue
                    yield self.sumexp[plate][well][org_id]
                else:               
                    for strain in sorted(self.sumexp[plate][well].keys()):  
                        yield self.sumexp[plate][well][strain]
                    
    def getAverageSinglePlates(self):
        '''
        Generator to the SinglePlates (average)
        '''
        for plate in self.sumexp:
            d = {}
            for well in self.sumexp[plate]:
                for strain in  self.sumexp[plate][well]:
                    if strain not in d:
                        d[strain] = SinglePlate()
                        d[strain].plate_id = plate
                        d[strain].strain = strain
                    d[strain].data[well] = self.sumexp[plate][well][strain]
            
            for strain in d:
                yield d[strain]
    
    def trim(self, trimTime = None):
        '''
        Set the maximum time for each well by using the lowest value in the
        experiment.
        Returns the trim time.

        If trimTime is set, that time will be used
        '''
        if trimTime is not None:
            mtime = trimTime
            if mtime > self.getMaxTime():
                logger.warning('Selected trim time > then max time (%f vs. %f)'%(
                                mtime, self.getMaxTime()))
        else:
            mtime = self.getMinTime()
        
        for w in self.getWells(False):
            to_del = set()
            for time in w.signals.keys():
                if time > mtime:
                    to_del.add(time)
            for time in to_del:
                del w.signals[time]
                
        return mtime
    
    def purgeReplicas(self, policy='keep-min', delta=1, replica=None):
        '''
        Analyze the replicas and remove the outliers using one of the policies
        
        keep-min --> keep the replicas around the minimum
        keep-max --> keep the replicas around the maximum
        for the above policies, the outliers are delta activity steps over or
        below the minimum/maximum activity
        
        keep-min-one --> keep the smaller replica
        keep-max-one --> keep the bigger replica
        
        replica --> remove a specific replica
        
        The mean activity value is then stored as a Well object inside sumexp
        
        The discarded wells are stored as biolog_ids in a set (discarded)
        '''
        # Check the provided policy
        if policy not in self.policies:
            logger.error('Policy not recognized %s'%policy)
            return False
        
        if policy == 'replica':
            if replica not in self.getDistinctReplica():
                logger.error('Replica %d not present'%replica)
                return False
            
            for plate in self.experiment:
                for well in self.experiment[plate]:
                    for strain in self.experiment[plate][well]:
                        if replica in self.experiment[plate][well][strain]:
                            self.discarded.add((plate, well,
                                                strain, replica))
                            
                            del self.experiment[plate][well][strain][replica]
                            
                            # TODO: simplify here
                            rem_p = filter(lambda x: x.replica == replica,
                                        self.plates[plate].strains[strain])[0]
                            del rem_p.data[well]
                            
                            logger.debug('Purged %s %s %s %d'%(plate, well,
                                                           strain, replica))
            
            return True
        
        for plate in self.experiment:
            for well in self.experiment[plate]:
                for strain in self.experiment[plate][well]:  
                    reps = self.experiment[plate][well][strain].values()
                    act = np.array([x.activity for x in reps])
                    
                    if policy == 'keep-min' or policy == 'keep-min-one':
                        m = act.min()
                    
                    if policy == 'keep-max' or policy == 'keep-max-one':
                        m = act.max()
                        
                    # If a keep-one policy is on, choose the replica that
                    # matches the policy as close as possible
                    # (i.e. keep-min-one --> replica with smaller average signal)    
                    if policy == 'keep-min-one' or policy == 'keep-max-one':
                        candidates = filter(lambda x: x.activity == m, reps)
                        if len(candidates) == 1:
                            self.sumexp[plate][well][strain] = candidates[0]
                        elif len(candidates) == 0:
                            logger.critical('This shouldn\'t be possible!')
                            return False
                        else:
                            # Keep the best replica according to the policy
                            candidates = sorted(candidates, key=lambda x: x.area)
                            if policy == 'keep-min-one':
                                self.sumexp[plate][well][strain] = candidates[0]
                                for x in candidates[1:]:
                                    candidates.remove(x)
                            else:
                                self.sumexp[plate][well][strain] = candidates[-1]
                                for x in candidates[:-1]:
                                    candidates.remove(x)
                                                                
                    # Keep those replica distant at max delta from the
                    # minimum-maximum
                    if policy == 'keep-min' or policy == 'keep-max':
                        if policy == 'keep-min':
                            candidates = filter(lambda x: x.activity <= m + delta,
                                                reps)
                            if len(candidates) == 1:
                                self.sumexp[plate][well][strain] = candidates[0]
                            elif len(candidates) == 0:
                                logger.critical('This shouldn\'t be possible!')
                                return False
                            else:
                                # Keep the average activity
                                self.sumexp[plate][well][strain].activity = np.array(
                                                             [x.activity
                                                              for x in candidates]
                                                                 ).mean()
                        else:
                            candidates = filter(lambda x: x.activity >= m - delta,
                                                reps)
                            if len(candidates) == 1:
                                self.sumexp[plate][well][strain] = candidates[0]
                            elif len(candidates) == 0:
                                logger.critical('This shouldn\'t be possible!')
                                return False
                            else:
                                # Keep the average activity
                                self.sumexp[plate][well][strain].activity = np.array(
                                                             [x.activity
                                                              for x in candidates]
                                                                 ).mean()
                    
                    # Remove the outliers
                    for w in reps:
                        if w not in candidates:
                            self.discarded.add((plate, well,
                                                strain, w.replica))
                            
                            del self.experiment[plate][well][strain][w.replica]
                            
                            # TODO: simplify here
                            rem_p = filter(lambda x: x.replica == w.replica,
                                        self.plates[plate].strains[strain])[0]
                            del rem_p.data[well]
                            
                            logger.debug('Purged %s %s %s %d'%(plate, well,
                                                           strain, w.replica))
            
        self.purged = True
        return True
    
    def getDistinctReplica(self):
        '''
        Returns a set with all the distinct replica values
        '''
        rep = set()
        for w in self.getWells(False):
            rep.add(w.replica)
        return rep
    
    def getDistinctActivity(self):
        '''
        Returns a set with all the distinct AV values
        '''
        av = set()
        for w in self.getWells(False):
            av.add(w.activity)
        return av
    
    def getMaxActivity(self):
        '''
        Get the maximum activity
        Which is also the number of clusters used...
        '''
        return max([w.activity for w in self.getWells(False)])
    
    def getMaxParam(self, param):
        '''
        Get the maximum value for a certain parameter
        '''
        return max([getattr(w, param) for w in self.getWells()])
    
    def getMinTime(self):
        '''
        Get the minimum time
        '''
        return min([w.getMaxTime() for w in self.getWells(False)])
    
    def getMaxTime(self):
        '''
        Get the maximum time
        '''
        return max([w.getMaxTime() for w in self.getWells(False)])
    
    def setMaxParams(self):
        '''
        Find the maximum value of each parameter
        '''
        w = Well('dummy', 'dummy')
        self.maxParams['zero'] = {}
        self.maxParams['nonzero'] = {}
        for param in w.params:
            self.maxParams['zero'][param] = self.maxParams.get(param, 0)
            self.maxParams['nonzero'][param] = self.maxParams.get(param, 0)
            
        for w in self.getWells(False):
            if self.zero and w.plate_id in self.zeroPlates:
                z = 'zero'
            else:
                z = 'nonzero'
                
            for param in w.params:
                if getattr(w, param) > self.maxParams[z][param]:
                    self.maxParams[z][param] = getattr(w, param)
    
    def normalizeParam(self, param, value, zero=False):
        '''
        Take a parameter and return its normalization
        '''
        if self.maxParams == {}:
            self.setMaxParams()
        
        if zero:
            z = 'zero'
        else:
            z = 'nonzero'
        
        try:  
            return float(value)/float(self.maxParams[z][param])
        except ZeroDivisionError:
            return 0
    
    def _prepareClusters(self):
        if self.zero:
            dWells = {'zero':[],
                      'nonzero':[]}
            dParams = {'zero':[],
                       'nonzero':[]}
        else:
            dWells = {'nonzero':[]}
            dParams = {'nonzero':[]}
        
        for param in self.getWells():
            if self.zero and param.plate_id in self.zeroPlates:
                dWells['zero'].append(param)
                dParams['zero'].append([self.normalizeParam('max', removeNegatives(purgeNAN(param.max)), True),
                            self.normalizeParam('area', removeNegatives(purgeNAN(param.area)), True), 
                            self.normalizeParam('height', removeNegatives(purgeNAN(param.height)), True),
                            self.normalizeParam('lag', removeNegatives(purgeNAN(param.lag)), True),
                            self.normalizeParam('slope', removeNegatives(purgeNAN(param.slope)), True)])
            else:
                dWells['nonzero'].append(param)
                dParams['nonzero'].append([self.normalizeParam('max', removeNegatives(purgeNAN(param.max))),
                           self.normalizeParam('area', removeNegatives(purgeNAN(param.area))),
                           self.normalizeParam('height', removeNegatives(purgeNAN(param.height))),
                           self.normalizeParam('lag', removeNegatives(purgeNAN(param.lag))),
                           self.normalizeParam('slope', removeNegatives(purgeNAN(param.slope)))])
        
        return dParams, dWells
    
    def elbowTest(self, nrange=range(2, 13)):
        '''
        Perform an elbow test on the k-means clustering
        nrange should be a list with each n going to be used in the clusterization
        '''
        from ductape.phenome.clustering import _kmeans
        from ductape.phenome.clustering import getSseKmeans
        from ductape.phenome.clustering import plotElbow
        
        params_labels = ['max', 'area', 'height', 'lag', 'slope']
        
        dParams, dWells = self._prepareClusters()
        if self.zero:
            params = dParams['zero'] + dParams['nonzero']
        else:
            params = dParams['nonzero']
        
        X = np.array( params )
        
        dSse = {}
        for n_clust in nrange:
            logger.info('K-means clusterization (k=%d)'%n_clust)
            k_means = _kmeans(X, n_clust)
            dSse[n_clust] = getSseKmeans(k_means, X)
        
        plotElbow(dSse, params_labels)
    
    def clusterize(self, save_fig=False, n_clusters=10):
        '''
        Perform the biolog data clusterizzation
        The data is divided in two chunks if Zero subtraction has been done
        '''
        from ductape.phenome.clustering import mean, kmeans, plotClusters
        
        params_labels = ['max', 'area', 'height', 'lag', 'slope']
        
        dParams, dWells = self._prepareClusters()
        
        # Add some fake wells with no signal to make sure we will got a 
        # "zero cluster"
        if self.zero  and len(dParams['zero']) >= 1:
            for i in range(1,97):
                who = self.well()
                who.replica = 0
                who.plate_id = 'fake'
                who.well_id = 'fake'
                who.strain = 'fake'
                dWells['zero'].append(who)
                dParams['zero'].append([0.0, 0.0, 0.0, 0.0, 0.0])
        
        if len(dParams['nonzero']) >= 1:
            for i in range(1,97):
                who = self.well()
                who.replica = 0
                who.plate_id = 'fake'
                who.well_id = 'fake'
                who.strain = 'fake'
                dWells['nonzero'].append(who)
                dParams['nonzero'].append([0.0, 0.0, 0.0, 0.0, 0.0])
        
        # Perform the actual clusterizzations
        
        # Zero subtracted signals
        
        # "Control" MeanShift
        # If we will get 1 cluster, we have a real "flat" experiment
        # Fixed KMeans to get an activity scale
        if self.zero and len(dParams['zero']) >= 1:
            xZero = [x for x in dParams['zero']]
            m_z_labels = mean(xZero)
            k_z_labels = kmeans(xZero, n_clusters)
        
        if self.zero  and len(dParams['zero']) >= 1:
            m_z_nclusters = len(np.unique(m_z_labels))
            k_z_nclusters = len(np.unique(k_z_labels))
            if m_z_nclusters == 1:
                logger.warning('The zero-subtracted subset seems to have no activity!')
                self.setNoActivity()
            else:
                # Order the clusters by average area
                mArea = []
                dP = np.array(dParams['zero'])
                for k in range(k_z_nclusters):
                    my_members = k_z_labels == k
                    mArea.append((k, dP[my_members, 1].mean()))
                mArea = sorted(mArea, key=lambda x: x[1])
                
                dConvert = {}
                i = 0
                for t in mArea:
                    dConvert[t[0]] = i
                    i += 1
                
                # List for intermediate plotting
                k_z_activity = []
                
                for i in range(len(k_z_labels)):
                    who = dWells['zero'][i]
                    who.activity = dConvert[k_z_labels[i]]
                    k_z_activity.append(dConvert[k_z_labels[i]])
                
                # Intermediate plot
                if save_fig:
                    plotClusters(xZero, k_z_activity,
                                 params=params_labels,
                                 method='kmeans', prefix='zero')
        
        # Non subtracted signals
        
        # "Control" MeanShift
        # If we will get 1 cluster, we have a real "flat" experiment
        # Fixed KMeans to get an activity scale
        
        if len(dParams['nonzero']) >= 1:
            xNonZero = [x for x in dParams['nonzero']]
            m_nz_labels = mean(xNonZero)
            k_nz_labels = kmeans(xNonZero, n_clusters)
        
        if len(dParams['nonzero']) >= 1:
            m_nz_nclusters = len(np.unique(m_nz_labels))
            k_nz_nclusters = len(np.unique(k_nz_labels))
            if m_nz_nclusters == 1:
                logger.warning('The nonzero-subtracted subset seems to have no activity!')
                self.setNoActivity()
            else:
                # Order the clusters by average area
                mArea = []
                dP = np.array(dParams['nonzero'])
                for k in range(k_nz_nclusters):
                    my_members = k_nz_labels == k
                    mArea.append((k, dP[my_members, 1].mean()))
                mArea = sorted(mArea, key=lambda x: x[1])
                
                dConvert = {}
                i = 0
                for t in mArea:
                    dConvert[t[0]] = i
                    i += 1
                
                # List for intermediate plotting
                k_nz_activity = []
                
                for i in range(len(k_nz_labels)):
                    who = dWells['nonzero'][i]
                    who.activity = dConvert[k_nz_labels[i]]
                    k_nz_activity.append(dConvert[k_nz_labels[i]])
                    
                # Intermediate plot
                if save_fig:
                    plotClusters(xNonZero, k_nz_activity,
                                 params=params_labels,
                                 method='kmeans', prefix='nonzero')
    
    def _saveFigure(self, fig, title='', name='', svg=False):
        '''
        Fix and save a multiaxes figure
        '''
        fig.suptitle(title, size='large')
        
        cNorm  = colors.Normalize(vmin=0, vmax=self.getMaxActivity())
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cm.RdYlGn)
        scalarMap.set_array(np.array(range(10)))
        cax = fig.add_axes([0.925, 0.2, 0.03, 0.6])
        cax.text(0.50, 1.01, 'Activity', size=10, ha='center')       
        plt.colorbar(scalarMap, cax=cax)
        
        if svg:
            ftype = 'svg'
        else:
            ftype = 'png'
        
        plt.savefig('%s.%s'%(name,ftype))
        
        logger.info('Saved "%s" graph (%s.%s)'%(title, name, ftype))
        
        plt.clf()
    
    def plot(self, svg=False):
        '''
        Go for the overall plots!
        Colored according to the activity.
        '''
        fig = plt.figure(figsize=(12,6))
        
        logger.debug('Plotting overall Zero wells')
        ax = fig.add_subplot(1,2,1)
        self._plot(self.getZeroWells(params=False), 'ZeroPlot',
                  'Zero subtracted wells', svg, ax)
        
        logger.debug('Plotting overall NoZero wells')
        ax = fig.add_subplot(1,2,2)
        self._plot(self.getNoZeroWells(params=False), 'NoZeroPlot',
                   'NoZero subtracted wells', svg, ax)
        
        self._saveFigure(fig, 'Overall plot', 'Overall', svg)
        plt.clf()
        
        fig = plt.figure(figsize=(24,12))
        axid = 1
        
        for categ, wells in self.getCategoryWells(params=False):
            ax = fig.add_subplot(2,4,axid)
            
            logger.debug('Plotting overall %s wells'%categ)
            self._plot(wells, 'CategPlot_%s'%categ,
                   '%s'%categ, svg, ax)
            
            axid += 1
        
        self._saveFigure(fig, 'Overall plot (categories)', 'OverallCateg', svg)
    
    def _plot(self, iterwells, name, description='Overall plot', svg=False,
              axis=None):
        '''
        Plot all the wells in a single plot!
        Coloured according to the activity. 
        '''
        from ductape.common.utils import rangeColors
        
        if not axis:
            # Figure creation
            fig = plt.figure(figsize=(8,8))
            ax = fig.add_subplot(111)
        else:
            ax = axis
            
        ax.set_xlabel('Hour', size='small')
        ax.set_ylabel('Signal', size='small')
        
        color = rangeColors(0, self.getMaxActivity(),
                             cm.RdYlGn(np.arange(0,256)))
        
        counter = 0
        maxsig = 0.0
        maxtime = 0.0
        for w in iterwells:
            counter += 1
            
            if not w.compressed:
                w.compress()
            if not w.smoothed:
                # More aggressive smooth
                try:
                    w.smooth(30)
                except:
                    w.smooth()
            
            times = sorted(w.signals.keys())
            ax.plot(times, [w.signals[t] for t in times], color=color[w.activity],
                        rasterized=True)
            
            msig = max(w.signals.values())
            if msig > maxsig:
                maxsig = msig  
               
            mtime = max(w.signals.keys())
            if mtime > maxtime:
                maxtime = mtime
        
        if counter == 0:
            return
        
        ax.set_ylim(0,maxsig)
        ax.set_xlim(0,maxtime)
        
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))
        
        ax.set_title('%s'%description, size='small')
        
        if not axis:
            cNorm  = colors.Normalize(vmin=0, vmax=self.getMaxActivity())
            scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cm.RdYlGn)
            scalarMap.set_array(np.array(range(10)))
            cax = fig.add_axes([0.925, 0.2, 0.03, 0.6])
            cax.text(0.50, 1.05, 'Activity', size=10, ha='center')       
            plt.colorbar(scalarMap, cax=cax)
            
            if svg:
                ftype = 'svg'
            else:
                ftype = 'png'
            
            plt.savefig('%s.%s'%(name,ftype))
            plt.clf()

class BiologParser(object):
    '''
    Abstract class for parsing of PM data files
    '''
    _start = 'Data File'
    _plate = 'Plate Type'
    _strainType = 'Strain Type'
    _sample = 'Sample Number'
    _strainName = 'Strain Name'
    _strainNumber = 'Strain Number'
    _other = 'Other'
    _dataStart = 'Hour'
    
    _platesPrefix = 'PM'
    
    def __init__(self, infile, validPlates=[]):
        # Biolog
        self.file = infile
        
        # Results
        self.plates = []
    
    def parse(self):
        try:
            self.parseOPM()
        except Exception as e:
            logger.warning('YAML/OPM parsing failed!')
            logger.debug('%s'%e)
            try:
                self.parseCSV()
            except Exception as e:
                logger.error('CSV parsing failed!')
                logger.debug('%s'%e)
                return False
        return True
    
    def parseCSV(self):
        plate = None
        data = False
        wells = []
        
        tblreader = csv.reader(open(self.file, 'rbU'), delimiter=',',
                               quotechar='"')
        for line in tblreader:
            if len(line) < 2:
                continue
            elif self._start in line[0].strip():
                # Do we have to save the old plate?
                if plate:
                    self.plates.append(plate)
                data = False
                wells = []
                plate = SinglePlate()
            elif self._plate in line[0].strip():
                plateID = line[1].strip()
                
                # Parse also non-standard plate IDs
                if not plateID.startswith(self._platesPrefix):
                    logger.warning('Non-standard plate ID found (%s)'%plateID)
                    logger.warning('Plate IDs should start with %s'%self._platesPrefix)
                    plate.plate_id = plateID
                    continue
                    
                # Simplify the plates IDs, removing letters, as opm does
                pID = plateID[2:]
                while len(pID) > 0:
                    try:
                        int(pID)
                        break
                    except ValueError:
                        pID = pID[:-1]
                
                # No luck
                if len(pID) == 0:
                    logger.warning('Non-standard plate ID found (%s)'%plateID)
                    plate.plate_id = plateID
                    continue
                elif int(pID) < 0:
                    logger.warning('Non-standard plate ID found (%s)'%plateID)
                    plateID = self._platesPrefix + abs(int(pID))
                    logger.warning('Going to use this ID (%s)'%plateID)
                    plate.plate_id = plateID
                    continue
                
                plateID = self._platesPrefix + '%02d'%int(pID)
                plate.plate_id = plateID
                    
            elif self._strainType in line[0].strip():
                if not plate:continue
                plate.strainType = line[1].strip()
            elif self._sample in line[0].strip():
                if not plate:continue
                plate.sample = line[1].strip()
            elif self._strainNumber in line[0].strip():
                if not plate:continue
                plate.strainNumber = line[1].strip()
            elif self._strainName in line[0].strip():
                if not plate:continue
                plate.strainName = line[1].strip()
            elif self._other in line[0].strip():
                if not plate:continue
                plate.other = line[1].strip()
            elif self._dataStart in line[0].strip():
                if not plate:continue
                data = True
                for i in range(len(line)):
                    if i == 0:continue
                    x = line[i]
                    if x == '':continue
                    plate.data[x.strip()] = Well(plate.plate_id, x.strip())
                    plate._idx[i] = x.strip()
                    wells.append(x.strip())
            elif data:
                if not plate:continue
                # Workaround for bad-formatted files
                try: float(line[0])
                except:
                    logger.debug('Could not parse this line from biolog file (%s)'%line)
                    continue
                #
                
                time = float(line[0])
                for i in range(len(line)):
                    if i == 0:continue
                    x = line[i]
                    if x == '':continue
                    well = plate._idx[i]
                    plate.data[well].addSignal(time, float(x))
        
        # The last plate should be saved as well!
        if plate and plate not in self.plates:
            self.plates.append(plate)
        
        return True
    
    def parseOPM(self):
        import yaml
        
        data = yaml.load(open(self.file))
        
        # We can have one single plate or several
        # we need to discriminate
        try:
            data.keys()
            data = [data]
        except:pass
            
        for pobj in data:
            plate = SinglePlate()
            
            # General plate attributes
            plateID = pobj['csv_data'][self._plate]
                
            # Parse also non-standard plate IDs
            if not plateID.startswith(self._platesPrefix):
                logger.warning('Non-standard plate ID found (%s)'%plateID)
                logger.warning('Plate IDs should start with %s'%self._platesPrefix)
                plate.plate_id = plateID
            else:
                    
                # Simplify the plates IDs, removing letters, as opm does
                pID = plateID[2:]
                while len(pID) > 0:
                    try:
                        int(pID)
                        break
                    except ValueError:
                        pID = pID[:-1]
                
                # No luck
                if len(pID) == 0:
                    logger.warning('Non-standard plate ID found (%s)'%plateID)
                    plate.plate_id = plateID
                elif int(pID) < 0:
                    logger.warning('Non-standard plate ID found (%s)'%plateID)
                    plateID = self._platesPrefix + abs(int(pID))
                    logger.warning('Going to use this ID (%s)'%plateID)
                    plate.plate_id = plateID
                else:             
                    plateID = self._platesPrefix + '%02d'%int(pID)
                    plate.plate_id = plateID
                
            plate.strainType = pobj['csv_data'][self._strainType]
            plate.sample = pobj['csv_data'][self._strainType]
            plate.strainNumber = pobj['csv_data'][self._strainNumber]
            plate.strainName = pobj['csv_data'][self._strainName]
            plate.other = pobj['csv_data'][self._other]
            
            # Well signals
            # we assume that they are always there
            times = pobj['measurements']['Hour']
            for wid in pobj['measurements']:
                if wid == 'Hour':continue
                
                plate.data[wid] = Well(plate.plate_id, wid)
                for i in range(len(times)):
                    plate.data[wid].addSignal(times[i],
                                              pobj['measurements'][wid][i])
            
            # Curve parameters
            # Do we have them?
            if 'aggregated' not in pobj:
                self.plates.append(plate)
                continue
            
            # Collect the software source 
            if 'aggr_settings' in pobj and 'software' in pobj['aggr_settings']:
                soft = pobj['aggr_settings']['software']
                for wid in plate.data:
                    plate.data[wid].source = soft
            
            # Collect the curve parameters
            for wid in pobj['aggregated']:
                plate.data[wid].max = nullifyNAN(pobj['aggregated'][wid]['A'])
                plate.data[wid].area = nullifyNAN(pobj['aggregated'][wid]['AUC'])
                plate.data[wid].lag = nullifyNAN(pobj['aggregated'][wid]['lambda'])
                plate.data[wid].slope = nullifyNAN(pobj['aggregated'][wid]['mu'])
            
            self.plates.append(plate)
            
        return True

class BiologZero(object):
    '''
    Class BiologZero
    Takes a list of SinglePlate objects and subtract the zero signal
    Two modes: 
        - normal (subtract the first well)
        - blank plate (subtract the zero plate --> zero list)
    plates [SinglePlate] --> well_id --> Well
    '''
    def __init__(self, data, blank=False, blankData=[], 
                 forceZero = True, zeroPlates=[], controlWells={},
                 zeroWells={}):
        # Biolog
        self.data = data
        self.blank = bool(blank)
        self.blankData = blankData
        self.forceZero = bool(forceZero)
        
        # User provided rules for control subtraction
        self.zeroPlates = set(zeroPlates)
        self.controlWells = controlWells
        self.zeroWells = zeroWells
        
        # Results
        self.plates = []
    
    def _zeroNormal(self, plate):
        '''
        Normal zero subtraction
        For some plates the first well is a the negative control
        '''
        if plate.plate_id in self.zeroPlates:
            if(plate.plate_id not in self.controlWells or
                    plate.plate_id not in self.zeroWells):
                logger.warning('Missing control well information: '+
                           'zero subtraction on plate %s aborted'%plate.plate_id)
                return
            
            for well in plate.data:
                # Is it a control well?
                if well in self.controlWells[plate.plate_id]:
                    continue
                # Get the specific control well for this well
                zero = plate.data[ self.zeroWells[plate.plate_id][well] ]
                # We assume that wells from the same plate
                # will end at the same time
                for hour in sorted(zero.signals.keys()):
                    plate.data[well].signals[hour] -= zero.signals[hour]
                    # Values below zero are forced to zero
                    if self.forceZero and plate.data[well].signals[hour] <= 0:
                        plate.data[well].signals[hour] = 0.1
                    
            # Last step: put the control wells to zero
            for zerowell in self.controlWells[plate.plate_id]:
                zero = plate.data[zerowell]
                for hour in zero.signals.keys():
                    zero.signals[hour]=0
                    
    def _zeroBlank(self, plate):
        '''
        Blank plate zero subtraction
        We assume that the user will provide us all the plates we need
        '''
        # Search fir the right plate
        found = False
        for zplate in self.blankData:
            if zplate.plate_id != plate.plate_id:
                continue
            found = True
            # Blank plates warning
            for well in plate.data:
                # We CANNOT assume that wells from the same plate
                # will end at the same time
                # in case of errors the missing values will be set to zero
                for hour in sorted(zplate.data[well].signals.keys()):
                    try:
                        plate.data[well].signals[hour] -= zplate.data[well].signals[hour]
                        # Values below zero are forced to zero
                        if self.forceZero and plate.data[well].signals[hour] <= 0:
                            plate.data[well].signals[hour] = 0.1
                    except:
                        logger.debug('Time %f present in blank plate was not'%(hour)+ 
                                        ' found on plate %s, signal was forced to'%(plate.plate_id)+
                                        ' zero')
                # Reset those hours that are not present in the blank plate
                for hour in sorted(plate.data[well].signals.keys()):
                    if hour not in zplate.data[well].signals:
                        logger.debug('Time %f present in plate %s was not'%(hour, plate.plate_id)+
                                        ' found on blank plate, signal was forced to'+
                                        ' zero')
                        plate.data[well].signals[hour] = 0.1
            
        if not found:
            logger.warning('Blank plate zero subtraction: could not find'+
                                ' blank plate %s'%plate.plate_id)
    
    def zeroSubTract(self):
        '''
        Zero subtraction
        '''
        if self.blank:
            logger.info('Zero subtraction with blank plate')
        else:
            logger.info('Normal zero subtraction')
        
        for plate in self.data:
            if self.blank:
                self._zeroBlank(plate)
            else:
                self._zeroNormal(plate)
                
            plate.zero = True
            for well in plate.getWells():
                well.zero = True
                
        self.plates = self.data
                
        return True
        
class BiologPlot(CommonThread):
    '''
    Class BiologPlot
    Takes a list of SinglePlate objects and creates some plots
    '''
    _statusDesc = {0:'Not started',
                1:'Making room',
                2:'Preparing data',
                3:'Preparing data (average)',
                4:'Preparing plots',
                5:'Creating the legend',
                6:'Creating plates plots',
                7:'Creating single plots',
                8:'Creating Heat maps'}
    
    _substatuses = [2,3,5,6,7,8]
    
    def __init__(self, data, avgdata = [],
                 expname = 'exp', 
                 plateNames = {}, wellNames = {},
                 colors = {}, smooth = True, window = 11,
                 compress = 0,
                 maxsig = None,
                 plotPlates=True, plotAll=False, plotActivity=True,
                 order = [], category = {},
                 svg=False,
                 plate=None,
                 well=None,
                 queue=Queue.Queue()):
        CommonThread.__init__(self,queue)
        # Biolog
        self.data = data
        self.avgdata = avgdata
        
        self.expname = expname
        
        # Plot parameters
        self.colors = colors
        self.plateNames = plateNames
        self.wellNames = wellNames
        self.smooth = bool(smooth)
        self.window = int(window)
        self.compress = int(compress)
        if maxsig:
            self.maxsig = float(maxsig)
        else:
            self.maxsig = None
        self.plotPlates = bool(plotPlates)
        self.plotAll = bool(plotAll)
        self.plotActivity = bool(plotActivity)
        self.order = order
        self.category = category
        self.svg = bool(svg)
        
        # Single plate/well
        self.splate = plate
        self.swell = well
        
        # Results
        # Plate_id --> Plate
        self.results = {}
        self.avgresults = {}
        # Single well plot
        self.well = None
    
    def makeRoom(self,location=''):
        '''
        Creates a tmp directory in the desired location
        '''
        try:
            path = os.path.abspath(location)
            path = os.path.join(path, 'tmp')
            try:os.mkdir(path)
            except:pass
            path = os.path.join(path, 'plots')
            try:os.mkdir(path)
            except:pass
            path = os.path.join(path, self.expname)
            self._room = path
            try:os.mkdir(path)
            except:pass
            
            # Prepare the categories directories
            if len(self.category) > 0:
                for c in set(self.category.values()):
                    path = os.path.join(self._room, c)
                    try:os.mkdir(path)
                    except:pass
        except:
            logger.debug('Temporary directory creation failed! %s'
                          %path)
    
    def getPlot(self, plate_id, well_id):
        '''
        A specific well is plotted and assigned to attribute well
        To save memory just one well can be plotted at any time
        Returns True/False
        '''
        if plate_id not in self.results:
            logger.warning('Plate %s was not found!'%plate_id)
            return False
        if well_id not in self.results[plate_id].wells:
            logger.warning('Well %s was not found!'%well_id)
            return False
        
        if self.well:
            self.well.clf()
            self.well = self.results[plate_id].plotWell(well_id, self.well)
        else:
            self.well = self.results[plate_id].plotWell(well_id)
        
        return True
    
    def run(self):
        if self.splate == None:
            self.updateStatus()
            self.makeRoom()
            self.startCleanUp()
            self.makeRoom()
        else:
            self.updateStatus(send=False)
        
        if self.killed:
            return
        
        if self.splate == None:
            self._maxsubstatus = len(self.data)
        else:
            self._maxsubstatus = len(filter(lambda x: x.plate_id == self.splate,
                                            self.data))
        self.updateStatus()
        for plate in self.data:
            
            # Single plate mode
            if self.splate != None and plate.plate_id != self.splate:
                continue
            #
                
            self._substatus += 1
            self.updateStatus(True)
            logger.debug('Adding plate %s'%plate.plate_id)
            
            if plate.plate_id in self.plateNames:
                plate_name = self.plateNames[plate.plate_id]
            else:
                plate_name = ''
            if plate.plate_id not in self.results:
                self.results[plate.plate_id] = Plate(plate.plate_id, 
                                                    plate_name = plate_name,
                                                    smooth = self.smooth,
                                                    window = self.window,
                                                    compress = self.compress,
                                                    maxsig = self.maxsig)
            
            self.results[plate.plate_id].addData(plate.strain, plate)
            # Sanity check
            if plate.strain not in self.colors:
                logger.error('Color code for strain %s is missing'%plate.strain)
                return
            #
            self.results[plate.plate_id].setColor(plate.strain,
                                                  self.colors[plate.strain])
        self.resetSubStatus()
        
        if self.killed:
            return
        
        if self.splate == None:
            self._maxsubstatus = len(self.avgdata)
        else:
            self._maxsubstatus = 1
        if self.swell == None:
            self.updateStatus()
        else:
            self.updateStatus(send=False)
        for plate in self.avgdata:
            
            # Single well mode
            if self.swell != None:
                break
            # Single plate mode
            if self.splate != None and plate.plate_id != self.splate:
                continue
            #
        
            self._substatus += 1
            self.updateStatus(True)
            logger.debug('Adding average plate %s'%plate.plate_id)
            
            if plate.plate_id in self.plateNames:
                plate_name = self.plateNames[plate.plate_id]
            else:
                plate_name = ''
            if plate.plate_id not in self.avgresults:
                self.avgresults[plate.plate_id] = Plate(plate.plate_id, 
                                                    plate_name = plate_name,
                                                    smooth = self.smooth,
                                                    window = self.window,
                                                    compress = self.compress,
                                                    maxsig = self.maxsig)
            
            self.avgresults[plate.plate_id].addData(plate.strain, plate)
            # Sanity check
            if plate.strain not in self.colors:
                logger.error('Color code for strain %s is missing'%plate.strain)
                return
            #
            self.avgresults[plate.plate_id].setColor(plate.strain,
                                                  self.colors[plate.strain])
        self.resetSubStatus()
        
        if self.killed:
            return
        
        self.updateStatus()
        for plate_id in sorted(self.results.keys()):
            logger.debug('Preparing plate %s'%plate.plate_id)
            self.results[plate_id].preparePlot()
        for plate_id in sorted(self.avgresults.keys()):
            logger.debug('Preparing average plate %s'%plate.plate_id)
            self.avgresults[plate_id].preparePlot()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Plot the legend
        if self.splate == None:
            self._maxsubstatus = len(self.results)
        else:
            self._maxsubstatus = 1
        self.updateStatus()
        for plate_id in sorted(self.results.keys()):
            self._substatus += 1
            self.updateStatus(True)
            logger.debug('Preparing legend for plate %s'%plate_id)
            
            self.results[plate_id].plotLegend(plate_id, strains=self.order)
            
            if self.svg:
                fformat = 'svg'
            else:
                fformat = 'png'
            
            if self.splate == None:
                if plate_id in self.category:
                    path = os.path.join(self._room,self.category[plate_id])
                else:
                    path = self._room
                fname = os.path.join(path,'%s_legend.%s'%(plate_id, fformat))
            else:
                # Single mode, save here
                fname = '%s_legend.%s'%(plate_id, fformat)
            
            if self.svg:
                self.results[plate_id].legend.savefig(fname, dpi=300)
            else:
                self.results[plate_id].legend.savefig(fname, dpi=150)
            self.results[plate_id].legend.clf()
            
            if self.killed:
                return
            
        self.resetSubStatus()
        
        if self.plotPlates:
            if self.splate == None:
                self._maxsubstatus = len(self.results)*96
            elif self.splate != None and self.swell == None:
                self._maxsubstatus = 96
            
            if self.splate != None and self.swell != None:
                self.updateStatus(send=False)
            else:
                self.updateStatus()
            for plate_id in sorted(self.results.keys()):
                if self.splate != None and self.swell != None:
                    break
                #
                
                logger.debug('Plotting plate %s'%plate_id)
                for i in self.results[plate_id].plotAll():
                    self._substatus += 1
                    self.updateStatus(True)
                    
                    if self.killed:
                        return
                
                if self.svg:
                    fformat = 'svg'
                else:
                    fformat = 'png'
                 
                if self.splate == None:  
                    # TODO: remember qgraphicspixmapitem for GUI clickable!
                    if plate_id in self.category:
                        path = os.path.join(self._room,self.category[plate_id])
                    else:
                        path = self._room
                    fname = os.path.join(path,'%s.%s'%(plate_id, fformat))
                else:
                    fname = '%s.%s'%(plate_id, fformat)
                
                if self.svg:
                    self.results[plate_id].figure.savefig(fname, dpi=300)
                else:
                    self.results[plate_id].figure.savefig(fname, dpi=150)
                self.results[plate_id].figure.clf()
        else:
            self.updateStatus(send=False)
        self.resetSubStatus()
        
        if self.killed:
            return
        
        if self.plotAll or (self.splate != None and self.swell != None):
            if self.swell == None:
                self._maxsubstatus = len(self.results)*96
            else:
                self._maxsubstatus = 1
            if self.splate != None and self.swell == None:
                self.updateStatus(send=False)
            else:
                self.updateStatus()
            for plate_id in sorted(self.results.keys()):
            
                # Single plate
                if self.splate != None and self.swell == None:
                    break
                #
            
                if plate_id in self.wellNames:
                    self.results[plate_id].addWellTitles(self.wellNames[plate_id])
                    
                for well_id in sorted(self.results[plate_id].wells):
                    # Single well
                    if self.swell != None and self.swell != well_id:
                        continue
                    #
                    
                    logger.debug('Plotting %s %s'%(plate_id, well_id))
                    
                    if not self.getPlot(plate_id, well_id):
                        self.sendFailure('Single plot creation failure')
                        return
                        
                    if self.svg:
                        fformat = 'svg'
                    else:
                        fformat = 'png'
                        
                    if self.swell == None:
                        if plate_id in self.category:
                            path = os.path.join(self._room,self.category[plate_id])
                        else:
                            path = self._room
                        fname = os.path.join(path,'%s_%s.%s'%(plate_id,
                                                         well_id, fformat))
                    else:
                        fname = '%s_%s.%s'%(plate_id, well_id, fformat)
                    
                    if self.svg:
                        self.well.savefig(fname, dpi=300)
                    else:
                        self.well.savefig(fname, dpi=150)
                    
                    self._substatus += 1
                    self.updateStatus(True)
                    
                    if self.killed:
                        return
        else:
            self.updateStatus(send=False)
        self.resetSubStatus()
        
        if self.killed:
            return
        
        if self.plotActivity:
            if self.splate == None:
                self._maxsubstatus = len(self.results)*96
            else:
                self._maxsubstatus = 96
            if self.swell == None:
                self.updateStatus()
            else:
                self.updateStatus(send=False)
            if len(self.avgresults) != 0:
                maxAct = max([p.getMaxActivity() for pid, p in self.avgresults.iteritems()])
            else:
                maxAct = 0
            for plate_id in sorted(self.avgresults.keys()):
                # Single well:
                if self.swell != None:
                    break
                #
                
                logger.debug('Plotting heatmap %s'%plate_id)
                
                for i in self.avgresults[plate_id].plotActivity(strains=self.order, maxAct=maxAct):
                    self._substatus += 1
                    self.updateStatus(True)
                    
                    if self.killed:
                        return
                    
                if self.svg:
                    fformat = 'svg'
                else:
                    fformat = 'png'
                
                if self.splate == None:   
                    if plate_id in self.category:
                        path = os.path.join(self._room,self.category[plate_id])
                    else:
                        path = self._room
                    fname = os.path.join(path,'%sheat.%s'%(plate_id, fformat))
                else:
                    fname = '%sheat.%s'%(plate_id, fformat)
                
                if self.svg:
                    self.avgresults[plate_id].heatfig.savefig(fname)
                else:
                    self.avgresults[plate_id].heatfig.savefig(fname, dpi=150)
                self.avgresults[plate_id].heatfig.clf()
        else:
            self.updateStatus(send=False)
        self.resetSubStatus()

class CalcParams(object):
    def __init__(self, well):
        self.well = well
    
    def __call__(self):
        if not self.well.max:
            try:
                logger.debug('Calculating parameters for %s - %s'%
                             (self.well.plate_id, self.well.well_id))
                if self.well.isParams():
                    self.well.calculateParams()
                else:
                    logger.debug('Parameters already present, ',
                                 'skipping parameters calculation')
            except:
                return False
        
        return True

class BiologCluster(CommonThread):
    '''
    Class BiologCluster
    '''
    _statusDesc = {0:'Not started',
               1:'Calculating parameters',
               2:'Clustering'}
    
    _substatuses = [1]
    
    def __init__(self,experiment,
                 save_fig_clusters=False, force_params=False, n_clusters=10,
                 elbow=False,
                 queue=Queue.Queue()):
        CommonThread.__init__(self,queue)
        # Experiment
        self.exp = experiment
        
        # Number of clusters?
        self.n_clusters = int(n_clusters)
        
        # Save clusters figure?
        self.save_fig = bool(save_fig_clusters)
        
        # Force parameters calculation (even if they are there already?)
        self.force = bool(force_params)
        
        # Elbow test instead of clusterization?
        self.elbow = bool(elbow)
        
    def calculateParams(self):
        wellcount = 0
        for w in self.exp.getWells(params=False):
            wellcount += 1
        
        self._maxsubstatus = wellcount
        
        for well in self.exp.getWells(params=False):
            logger.debug('Calculating parameters for %s - %s'%
                             (well.plate_id, well.well_id))
            self._substatus += 1
            self.updateStatus(sub=True)
            
            if not well.isParams() or self.force:
                well.calculateParams()
            else:
                logger.debug('Parameters already present, '
                                 'skipping parameters calculation')
        
        return True
    
    def run(self):
        self.updateStatus()
        if not self.calculateParams():
            self.sendFailure('Calculate parameters failure!')
            return
        self.resetSubStatus()
        
        if self.killed or self.elbow:
            return
        
        self.updateStatus()
        self.exp.clusterize(self.save_fig, self.n_clusters)

def getSinglePlates(binput, nonmean=False):
    '''
    Takes signals or wells from the storage and transforms them into SinglePlates
    NB it is a generator
    '''
    if len(binput) == 0:
        return
    
    if hasattr(binput[0], "times"):
        for splate in getSinglePlatesFromSignals(binput):
            yield splate
    else:
        for splate in getSinglePlatesFromParameters(binput, nonmean):
            yield splate
            
def getSinglePlatesFromSignals(signals):
    '''
    Takes a bunch of signals taken from the DB and returns a series of 
    SinglePlates objects
    NB it is a generator
    '''
    dExp = {}
    
    for well in signals:
        plate_id, well_id, org_id, replica = (well.plate_id, well.well_id,
                                              well.org_id, well.replica)
        
        lT = well.times.split('_')
        lS = well.signals.split('_')
        
        if plate_id not in dExp:
            dExp[plate_id] = {}
        if org_id not in dExp[plate_id]:
            dExp[plate_id][org_id] = {}
        if replica not in dExp[plate_id][org_id]:
            dExp[plate_id][org_id][replica] = SinglePlate()
            dExp[plate_id][org_id][replica].plate_id = plate_id
            dExp[plate_id][org_id][replica].strain = org_id
            dExp[plate_id][org_id][replica].replica = replica
        if well_id not in dExp[plate_id][org_id][replica].data:
            dExp[plate_id][org_id][replica].data[well_id] = Well(plate_id,
                                                                 well_id)
        
        for i in range(len(lT)):
            dExp[plate_id][org_id][replica].data[well_id].addSignal(float(lT[i]),
                                                                    float(lS[i]))
            
        # Add the activity - if present
        if hasattr(well, "activity"):
            dExp[plate_id][org_id][replica].data[well_id].activity = well.activity
        #
        
        # And the other parameters as well
        for param in Well('fake', 'fake').params + Well('fake', 'fake').otherparams:
            if hasattr(well, param):
                setattr(dExp[plate_id][org_id][replica].data[well_id],
                        param,
                        getattr(well, param, None))
        
    # Return all the SinglePlates objects 
    for orgs in dExp.itervalues():
        for replicas in orgs.itervalues():
            for splate in replicas.itervalues():
                yield splate

def getSinglePlatesFromParameters(wells, nonmean=False):
    '''
    Takes a bunch of wells taken from the DB and returns a series of 
    SinglePlates objects
    NB it is a generator
    '''
    dExp = {}
    
    if nonmean:
        for well in wells:
            plate_id = well.plate_id
            well_id = well.well_id
            org_id = well.org_id
            replica = well.replica
            
            if plate_id not in dExp:
                dExp[plate_id] = {}
            if org_id not in dExp[plate_id]:
                dExp[plate_id][org_id] = {}
            if replica not in dExp[plate_id][org_id]:
                dExp[plate_id][org_id][replica] = SinglePlate()
                dExp[plate_id][org_id][replica].plate_id = plate_id
                dExp[plate_id][org_id][replica].strain = org_id
                dExp[plate_id][org_id][replica].replica = replica
            if well_id not in dExp[plate_id][org_id][replica].data:
                dExp[plate_id][org_id][replica].data[well_id] = Well(plate_id,
                                                                     well_id)
            
            dExp[plate_id][org_id][replica].data[well_id].activity = well.activity
            
            for param in Well('fake', 'fake').params + Well('fake', 'fake').otherparams:
                setattr(dExp[plate_id][org_id][replica].data[well_id],
                        param,
                        getattr(well, param, None)) 
            
        # Return all the SinglePlates objects 
        for orgs in dExp.itervalues():
            for replicas in orgs.itervalues():
                for splate in replicas.itervalues():
                    yield splate
    else:
        for well in wells:
            plate_id = well.plate_id
            well_id = well.well_id
            org_id = well.org_id
            
            if plate_id not in dExp:
                dExp[plate_id] = {}
            if org_id not in dExp[plate_id]:
                dExp[plate_id][org_id] = {}
            if well_id not in dExp[plate_id][org_id]:
                dExp[plate_id][org_id][well_id] = []
            
            dExp[plate_id][org_id][well_id].append(well.activity)
        
        # Return all the SinglePlates objects
        # After the calculation of the mean activity index
        for plate_id in dExp:
            for org_id in dExp[plate_id]:
                splate = SinglePlate()
                splate.plate_id = plate_id
                splate.strain = org_id
                splate.replica = 0
                for well_id in dExp[plate_id][org_id]:
                    splate.data[well_id] = Well(plate_id, well_id)
                    splate.data[well_id].activity = np.array(dExp[plate_id][org_id][well_id]).mean()
                yield splate
                
def getPlates(signals, nonmean=False):
    '''
    Takes a bunch of signals taken from the DB and returns a series of 
    Plates objects
    NB it is a generator
    '''
    dExp = {}
    for splate in getSinglePlates(signals, nonmean):
        if splate.plate_id not in dExp:
            dExp[splate.plate_id] = Plate(splate.plate_id)
        dExp[splate.plate_id].addData(splate.strain, splate)
    
    for plate in dExp.itervalues():
        yield plate

def toOPM(p):
    d={}
    
    d['csv_data'] = {}
    d['csv_data']['Data file'] = ''
    d['csv_data']['File'] = ''
    d['csv_data']['Other'] = ''
    d['csv_data']['Plate Type'] = p.plate_id
    d['csv_data']['Position'] = ''
    d['csv_data']['Sample Number'] = p.strain
    d['csv_data']['Setup Time'] = ''
    d['csv_data']['Strain Name'] = p.strain
    d['csv_data']['Strain Number'] = p.strain
    d['csv_data']['Strain Type'] = ''
    
    d['metadata'] = []
    
    d['measurements'] = {}
    d['measurements']['Hour'] = []
    times = set()
    for wid in p.data:
        d['measurements'][wid] = []
        for hour in p.data[wid].signals:
            times.add(hour)
    
    for hour in sorted(times):
        d['measurements']['Hour'].append(hour)
        for wid in p.data:
            if hour in p.data[wid].signals:
                d['measurements'][wid].append(p.data[wid].signals[hour])
            # This shouldn't happen
            else:
                d['measurements'][wid].append(float('nan'))
      
    # Do we have some parameters?
    isaggr = set([p.data[x].isParams() for x in p.data])
    if True not in isaggr:
        return d        
                
    d['aggr_settings'] = {}
    
    # We may have more than one parameters sources!
    sources = set([p.data[x].source for x in p.data])
    if len(sources) > 1:
        logger.error('Cannot export plate %s (replica %s, %s), having %d'
                     %(p.plate_id, p.replica, p.strain, len(sources))+
                     ' different parameters sources')
        for s in sources:
            logger.error('Method: %s'%s)
        raise ValueError('Too many different parameters sources')
    source = sources.pop()
    
    d['aggr_settings']['software'] = source
    d['aggr_settings']['method'] = source
    d['aggr_settings']['options'] = {'Dummy':'Dummy'}
    
    d['aggregated'] = {}
    for wid in p.data:
        d['aggregated'][wid] = {}
        
        if p.data[wid].slope is not None:
            d['aggregated'][wid]['mu'] = p.data[wid].slope
        else:
            d['aggregated'][wid]['mu'] = '.na.real'
            
        if p.data[wid].lag is not None:
            d['aggregated'][wid]['lambda'] = p.data[wid].lag
        else:
            d['aggregated'][wid]['lambda'] = '.na.real'
        
        if p.data[wid].max is not None:
            d['aggregated'][wid]['A'] = p.data[wid].max
        else:
            d['aggregated'][wid]['A'] = '.na.real'
            
        if p.data[wid].area is not None:
            d['aggregated'][wid]['AUC'] = p.data[wid].area
        else:
            d['aggregated'][wid]['AUC'] = '.na.real'
        
        d['aggregated'][wid]['mu CI95 low'] = '.na.real'
        d['aggregated'][wid]['lambda CI95 low'] = '.na.real'
        d['aggregated'][wid]['A CI95 low'] = '.na.real'
        d['aggregated'][wid]['AUC CI95 low'] = '.na.real'
        d['aggregated'][wid]['mu CI95 high'] = '.na.real'
        d['aggregated'][wid]['lambda CI95 high'] = '.na.real'
        d['aggregated'][wid]['A CI95 high'] = '.na.real'
        d['aggregated'][wid]['AUC CI95 high'] = '.na.real'
        
    return d

def toYAML(plate):
    '''
    Take a SimplePlate object and return YAML strings 
    '''
    import yaml
    return yaml.safe_dump(toOPM(plate), default_flow_style=False)

def toJSON(plate):
    '''
    Take a SimplePlate object and return JSON strings 
    '''
    import json
    return json.dumps(toOPM(plate))

def nullifyNAN(value):
    '''
    Takes a value, if it is not convertible to float returns None
    '''
    try:
        return float(value)
    except:
        return None
   
def purgeNAN(value):
    if value is None:
        return 0
    
    import math
    if math.isnan( float(value) ):
        return 0
    
    return value

def removeNegatives(value):
    if value < 0.0:
        return 0
    
    return value
