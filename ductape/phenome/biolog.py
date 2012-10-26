#!/usr/bin/env python
"""
Biolog

Phenome library

Classes to handle Biolog data
"""
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
# Plate ID conversion and sanity checks

dPlates = {'PM 1-':'PM01', 'PM 2-A':'PM02A',
            'PM 3-B':'PM03B', 'PM 4-A':'PM04A',
            'PM 5-':'PM05', 'PM 6-':'PM06',
            'PM 7-':'PM07', 'PM 8-':'PM08',
            'PM 9-':'PM09', 'PM10-':'PM10',
            'PM11-C':'PM11C', 'PM12-B':'PM12B',
            'PM13-B':'PM13B', 'PM14-A':'PM14A',
            'PM15-B':'PM15B', 'PM16-A':'PM16A',
            'PM17-A':'PM17A', 'PM18-C':'PM18C',
            'PM19-':'PM19', 'PM20-B':'PM20B',
            'PM01-':'PM01', 'PM02-A':'PM02',
            'PM03-B':'PM03B','PM04-A':'PM04A',
            'PM05-':'PM05', 'PM06-':'PM06',
            'PM07-':'PM07','PM08-':'PM08',
            'PM09-':'PM09',
            # TODO: These are just guesses for now
            'PM21-D':'PM21D', 'PM22-C':'PM22C',
            'PM23-A':'PM23A', 'PM24-B':'PM24B'}

acceptedPlates = dPlates.values()

zeroPlates = ['PM01','PM02A','PM03B','PM04A','PM05',
              'PM06','PM07','PM08','PM03B','PM04A']
zeroWell = 'A01'

wellChars = 'ABCDEFGH'

def getPlatesOrder():
    return sorted(set(dPlates.values()))

def getOrder(plates=None):
    '''
    Generator of plate/well IDs
    If plates is provided as a list of plates IDs, only those plates are used
    '''
    if not plates:
        plates = getPlatesOrder()
        
    for pid in plates:
        for wh in wellChars:
            for wi in range(1,13):
                if wi < 10:
                    wid = '%s0%d'%(wh, wi)
                else:
                    wid = '%s%d'%(wh, wi)
                yield (pid, wid)

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
                       'plateu', 'slope', 'lag',
                       'area', 'v', 'y0']
        
        # Fitting model used
        self.model = None
        
        # Relative activity index
        self.activity = None
        
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
                           'Slope', 'Lag', 'Area'] )
        
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
                                            self.area]] )

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
    
    def smooth(self, window_len = 11, window_type = 'hanning',
               forceZero = True):
        '''
        Apply a smoothing algorithm to the signals
        Really useful for clearer plots and other features
        Available windows: 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
        '''
        if not self.smoothed:
            signals = [self.signals[hour] for hour in sorted(self.signals.keys())]
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
        
        if not self.compressed and not noCompress:
            self.compress()
        if not self.smoothed and not noSmooth:
            self.smooth(window_len=len(self.signals)/3, window_type='blackman')
            
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
        
    def getWell(self):
        '''
        Generator to the single wells
        '''
        for well_id, well in self.data.iteritems():
            yield well
        
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
            if not well.max:
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
            acts = acts.reshape(square, square)
        else:
            # Oh snap!
            cols = int(square)
            rows = len(acts) / float(cols)
            if rows.is_integer():
                acts = acts.reshape(rows, cols)
            else:
                # Some fake signals will be added
                rows = int(rows) + 1
                new = (cols * rows) - len(acts)
                acts = np.array( acts.tolist() + [np.nan for i in range(new)] )
                acts = acts.reshape(rows, cols)
        
        return acts
    
    def plotActivity(self, strains=[]):
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
                raise ValueError('Unknown strain(s) were provided (%s)'%
                                 ' '.join(unknown))
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
                    
            # TODO: here the vmax value is hard-coded
            ax.matshow(acts, cmap=cm.RdYlGn, vmin=0, vmax=9)
            
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
                raise ValueError('Unknown strain(s) were provided (%s)'%
                                 ' '.join(unknown))
        else:
            strains = sorted(self.strains.keys())
        
        array = np.array([strain for strain in strains])
        array = self.arrayReshape(array, strains)
        
        self.legend = plt.figure()
        ax = self.legend.add_subplot(111)
        
        # Show the array
        acts = np.array([np.nan for strain in strains])
        acts = self.arrayReshape(acts, strains)
        
        # Setup the legend
        ax.matshow(acts)
        for i in range(acts.shape[0] - 1):
            ax.axhline(y=i+0.5, color='black')
        
        for i in range(acts.shape[0] - 1):
            ax.axvline(x=i+0.5, color='black')
        
        # Plot the strains names in the proper order
        for i in range(len(array)):
            for j in range(len(array[i])):
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
        
        data.replica = len(self.strains[strain])
        
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
                 category = {}, categorder = []):
        self.exp_id = exp_id
        self.name = name
        
        self.zero = zero
        self.category = category
        self.categorder = categorder
        
        self.plates = {}
        for plate in plates:
            if not self._addPlate(plate):
                self.plates = {}
                break
            
        self.experiment = {}
        self.sumexp = {}
        self._organize()
        
        # Allowed policies for purging of replicas
        self.policies = ['keep-min', 'keep-max',
                         'keep-min-one', 'keep-max-one']
        
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
                    
                    act = []
                    for r in reps:
                        if self.experiment[pid][wid][org][r].activity is None:
                            break
                        act.append(self.experiment[pid][wid][org][r].activity)
                        
                    if len(act) > 0:
                        mean = np.array(act).mean()
                        self.sumexp[pid][wid][org].activity = mean
    
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
    
    def getZeroWells(self, params=True):
        '''
        Generator to get the zero-subtracted single wells
        if params is set to False, it just gives you the wells,
        otherwise it calculates them
        '''
        for well in self.getWells(params):
            if well.plate_id in zeroPlates:
                yield well
                
    def getNoZeroWells(self, params=True):
        '''
        Generator to get the nonzero-subtracted single wells
        if params is set to False, it just gives you the wells,
        otherwise it calculates them
        '''
        for well in self.getWells(params):
            if well.plate_id not in zeroPlates:
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
                if not well.max and params:
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
                    for well in plate.data:
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
    
    def purgeReplicas(self, policy='keep-min', delta=1):
        '''
        Analyze the replicas and remove the outliers using one of the policies
        
        keep-min --> keep the replicas around the minimum
        keep-max --> keep the replicas around the maximum
        for the above policies, the outliers are delta activity steps over or
        below the minimum/maximum activity
        
        keep-min-one --> keep the smaller replica
        keep-max-one --> keep the bigger replica
        
        The mean activity value is then stored as a Well object inside sumexp
        
        The discarded wells are stored as biolog_ids in a set (discarded)
        '''
        # Check the provided policy
        if policy not in self.policies:
            logger.error('Policy not recognized %s'%policy)
            return False
        
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
    
    def clusterize(self, save_fig=False):
        '''
        Perform the biolog data clusterizzation
        The data is divided in two chunks if Zero subtraction has been done
        '''
        from ductape.phenome.clustering import mean, kmeans, plotClusters
        
        if self.zero:
            dWells = {'zero':[],
                      'nonzero':[]}
            dParams = {'zero':[],
                       'nonzero':[]}
        else:
            dWells = {'nonzero':[]}
            dParams = {'nonzero':[]}
        
        params_labels = ['max', 'area', 'height', 'lag', 'slope']
        for param in self.getWells():
            if self.zero and param.plate_id in zeroPlates:
                dWells['zero'].append(param)
                dParams['zero'].append([param.max, param.area, 
                                        param.height, param.lag, param.slope])
            else:
                dWells['nonzero'].append(param)
                dParams['nonzero'].append([param.max, param.area,
                                           param.height, param.lag, param.slope])
        
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
            k_z_labels = kmeans(xZero)
        
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
            k_nz_labels = kmeans(xNonZero)
        
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
        
        cNorm  = colors.Normalize(vmin=0, vmax=9)
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
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = axis
            
        ax.set_xlabel('Hour', size='small')
        ax.set_ylabel('Signal', size='small')
        
        colors = rangeColors(0, 9, cm.RdYlGn(np.arange(0,256)))
        
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
            ax.plot(times, [w.signals[t] for t in times], color=colors[w.activity],
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
            if svg:
                ftype = 'svg'
            else:
                ftype = 'png'
                
            plt.savefig('%s.%s'%(name,ftype))
            plt.clf()

class BiologParser(object):
    '''
    Class BiologParser
    Takes a csv output file and returns a series of objects
    plates [SinglePlate] --> well_id --> Well
    There could be more than one organism for each biolog file!
    '''
    _start = 'Data File'
    _plate = 'Plate Type'
    _strainType = 'Strain Type'
    _sample = 'Sample Number'
    _strainName = 'Strain Name'
    _strainNumber = 'Strain Number'
    _other = 'Other'
    _dataStart = 'Hour'
    
    def __init__(self,infile):
        # Biolog
        self.file = infile
        # Results
        self.plates = []
        
    def parse(self):
        plate = None
        data = False
        wells = []
        
        tblreader = csv.reader(open(self.file, 'rbU'), delimiter=',')
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
                if line[1].strip() in acceptedPlates:
                    plate.plate_id = line[1].strip()
                else:
                    plate.plate_id = dPlates[line[1].strip()]
            elif self._strainType in line[0].strip():
                plate.strainType = line[1].strip()
            elif self._sample in line[0].strip():
                plate.sample = line[1].strip()
            elif self._strainNumber in line[0].strip():
                plate.strainNumber = line[1].strip()
            elif self._strainName in line[0].strip():
                plate.strainName = line[1].strip()
            elif self._other in line[0].strip():
                plate.other = line[1].strip()
            elif self._dataStart in line[0].strip():
                data = True
                for i in range(len(line)):
                    if i == 0:continue
                    x = line[i]
                    if x == '':continue
                    plate.data[x.strip()] = Well(plate.plate_id, x.strip())
                    plate._idx[i] = x.strip()
                    wells.append(x.strip())
            elif data:
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
                 forceZero = True):
        # Biolog
        self.data = data
        self.blank = bool(blank)
        self.blankData = blankData
        self.forceZero = bool(forceZero)
        # Results
        self.plates = []
    
    def _zeroNormal(self, plate):
        '''
        Normal zero subtraction
        For some plates the first well is a the negative control
        '''
        if plate.plate_id in zeroPlates:
            zero = plate.data[zeroWell]
            for well in plate.data:
                if well == zero.well_id:
                    continue
                # We assume that wells from the same plate
                # will end at the same time
                for hour in sorted(zero.signals.keys()):
                    plate.data[well].signals[hour] -= zero.signals[hour]
                    # Values below zero are forced to zero
                    if self.forceZero and plate.data[well].signals[hour] <= 0:
                        plate.data[well].signals[hour] = 0.1
                    
            # Last step: put the zero well to zero
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
            for well in plate.getWell():
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
        self.updateStatus()
        self.makeRoom()
        self.startCleanUp()
        self.makeRoom()
        
        if self.killed:
            return
        
        self._maxsubstatus = len(self.data)
        self.updateStatus()
        for plate in self.data:
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
        
        self._maxsubstatus = len(self.avgdata)
        self.updateStatus()
        for plate in self.avgdata:
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
        self._maxsubstatus = len(self.results)
        self.updateStatus()
        for plate_id in sorted(self.results.keys()):
            self._substatus += 1
            self.updateStatus(True)
            logger.debug('Preparing legend for plate %s'%plate_id)
            
            self.results[plate_id].plotLegend(plate_id, strains=self.order)
            
            if plate_id in self.category:
                path = os.path.join(self._room,self.category[plate_id])
            else:
                path = self._room
            fname = os.path.join(path,'%s_legend.png'%plate_id)
            self.results[plate_id].legend.savefig(fname, dpi=150)
            self.results[plate_id].legend.clf()
            
            if self.killed:
                return
            
        self.resetSubStatus()
        
        if self.plotPlates:
            self._maxsubstatus = len(self.results)*96
            self.updateStatus()
            for plate_id in sorted(self.results.keys()):
                logger.debug('Plotting plate %s'%plate_id)
                for i in self.results[plate_id].plotAll():
                    self._substatus += 1
                    self.updateStatus(True)
                    
                    if self.killed:
                        return
                    
                # TODO: remember qgraphicspixmapitem for GUI clickable!
                if plate_id in self.category:
                    path = os.path.join(self._room,self.category[plate_id])
                else:
                    path = self._room
                fname = os.path.join(path,'%s.png'%plate_id)
                self.results[plate_id].figure.savefig(fname, dpi=150)
                self.results[plate_id].figure.clf()
        else:
            self.updateStatus(send=False)
        self.resetSubStatus()
        
        if self.killed:
            return
        
        if self.plotAll:
            self._maxsubstatus = len(self.results)*96
            self.updateStatus()
            for plate_id in sorted(self.results.keys()):
                if plate_id in self.wellNames:
                    self.results[plate_id].addWellTitles(self.wellNames[plate_id])
                    
                for well_id in sorted(self.results[plate_id].wells):
                    logger.debug('Plotting %s %s'%(plate_id, well_id))
                    
                    if not self.getPlot(plate_id, well_id):
                        self.sendFailure('Single plot creation failure')
                        return
                    if plate_id in self.category:
                        path = os.path.join(self._room,self.category[plate_id])
                    else:
                        path = self._room
                    fname = os.path.join(path,'%s_%s.png'%(plate_id,
                                                                 well_id))
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
            self._maxsubstatus = len(self.results)*96
            self.updateStatus()
            for plate_id in sorted(self.avgresults.keys()):
                logger.debug('Plotting heatmap %s'%plate_id)
                
                for i in self.avgresults[plate_id].plotActivity(strains=self.order):
                    self._substatus += 1
                    self.updateStatus(True)
                    
                    if self.killed:
                        return
                if plate_id in self.category:
                    path = os.path.join(self._room,self.category[plate_id])
                else:
                    path = self._room
                fname = os.path.join(path,'%sheat.png'%plate_id)
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
                self.well.calculateParams()
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
                 save_fig_clusters=False,queue=Queue.Queue()):
        CommonThread.__init__(self,queue)
        # Experiment
        self.exp = experiment
        
        # Save clusters figure?
        self.save_fig = bool(save_fig_clusters)
    
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
            
            well.calculateParams()
        
        return True
    
    def run(self):
        self.updateStatus()
        if not self.calculateParams():
            self.sendFailure('Calculate parameters failure!')
            return
        self.resetSubStatus()
        
        if self.killed:
            return
        
        self.updateStatus()
        self.exp.clusterize(self.save_fig)

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
            
            for param in Well('fake', 'fake').params:
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
