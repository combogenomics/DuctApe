#!/usr/bin/env python
"""
Biolog

Pheenome library

Classes to handle Biolog data
"""
# TODO: this part must be handled somewhere else
import matplotlib
matplotlib.use('Agg')
#
from DuctApe.Common.CommonThread import CommonThread
from DuctApe.Phenome.fitting import fitData
from DuctApe.Common.utils import smooth, compress
import Queue
import csv
import logging
import matplotlib.pyplot as plt
import numpy as np
import os

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('Biolog')

################################################################################
# Plate ID conversion

dPlates = {'PM 1-':'PM01', 'PM 2-A':'PM02',
            'PM 3-B':'PM03B', 'PM 4-A':'PM04A',
            'PM 5-':'PM05', 'PM 6-':'PM06',
            'PM 7-':'PM07', 'PM 8-':'PM08',
            'PM 9-':'PM09', 'PM10-':'PM10',
            'PM11-C':'PM11C', 'PM12-B':'PM12B',
            'PM13-B':'PM13B', 'PM14-A':'PM14',
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

zeroPlates = ['PM01','PM02','PM03','PM04','PM05',
              'PM06','PM07','PM08','PM03B','PM04A']
zeroWell = 'A01'

################################################################################
# Classes

class BiologExp(object):
    '''
    Class BiologExp
    General data about a particular well to be stored
    '''
    def __init__(self, plate_id, well_id, org_id, replica, active = False,
               prediction = None, zero = False):
        self.plate_id = plate_id
        self.well_id = well_id
        self.org_id = org_id
        self.replica = int(replica)
        self.active = bool(active)
        self.prediction = prediction
        self.zero = bool(zero)

class PlateCarrier(object):
    '''
    Class PlateCarrier
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
        
        # Raw data --> well_id -> BiologRaw objects
        self.data = {}
        # Used internaly, index in Hours row --> well_id
        self._idx = {}
        
    def getMax(self):
        '''
        Maximum signal for the entire plate
        '''
        return max( [self.data[well].getMax() for well in self.data] )

class PlotCarrier(object):
    '''
    Class PlotCarrier
    Contains all the distinct strains data for a particular plate
    There can be more than one replica for each strain
    An averaged plate can be added
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
        
        # Figures
        self.figure = None
        
        self._figidx = 1
        
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
                logging.error('Color code for strain %s is missing!'%strain)
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
            try:
                for plate in plates:
                    strain_signals[strain].append( 
                            [plate.data[well_id].signals[hour]
                             for hour in self.times])
            except:
                logging.debug('Something missing: %s, %s, %f'%(
                                    strain, well_id, hour))
                pass
            
        return strain_signals
    
    def plotAll(self):
        # Preparatory steps
        if not self.times and not self.wells:
            self.preparePlot()
        
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
            if abs(96 % self._figidx - 12) <= 12:
                ax.set_xlabel(str(abs(96 % self._figidx - 12)))
            self._figidx += 1
            
            yield well_id, self.wells.index(well_id)
        
        self.fixFigure()
    
    def plotWell(self, well_id):
        '''
        Generates and returns a single well as a figure
        '''
        # Preparatory steps
        if not self.times and not self.wells:
            self.preparePlot()
        
        strain_signals = self._prepareSignals(well_id)
        
        # Figure creation
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        self._plot(well_id, strain_signals, ax)
        
        if well_id in self.wellNames:
            ax.set_title(' '.join( [well_id, self.wellNames[well_id]]))
        ax.set_xlabel('Hour')
        ax.set_ylabel('Signal')
        
        return fig
        
    def setColor(self, strain, color):
        '''
        Set the color code of a strain
        '''
        self.colors[strain] = color
    
    def addData(self, strain, data):
        '''
        Add a PlateCarrier object regarding a particular strain
        A check on the plate_id is performed!
        '''
        if data.plate_id != self.plate_id:
            logging.error('Expecting %s, got %s'%(self.plate_id,data.plate_id))
            return False
        
        if strain not in self.strains:
            self.strains[strain] = []
        self.strains[strain].append(data)
        
        return True

class BiologRaw(object):
    '''
    Class BiologRaw
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
        
    def addSignal(self,time,signal):
        self.signals[time] = signal
        
    def fillMissing(self, times):
        '''
        Given a times list, fills the missing values with 0
        '''
        for hour in sorted(times):
            if hour not in self.signals:
                self.signals[hour] = 0.1
    
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
            logging.warning('Plate %s, Well %s was already smoothed'%
                          (self.plate_id, self.well_id))
            
    def compress(self, span = 3):
        '''
        Reduce the amount of signals
        This function should be called BEFORE smooth
        '''
        if self.smoothed:
            logging.warning('Plate %s, Well %s should be smoothed AFTER compression'%
                          (self.plate_id, self.well_id))
        
        if not self.compressed:
            times = compress( sorted(self.signals.keys()), span = span)
            toremove = [t for t in self.signals.keys() if t not in times]
            for t in toremove:
                del self.signals[t]
            
            self.compressed = True
        else:
            logging.warning('Plate %s, Well %s was already compressed'%
                          (self.plate_id, self.well_id))
            
    def calculateParameters(self,
                            noCompress = False, noSmooth = False):
        '''
        Populates the parameters values for the experiment
        By default compression and smoothing are applied to save some time
        '''
        if not self.compressed and not noCompress:
            self.compress()
        if not self.smoothed and not noSmooth:
            self.smooth(window_len=41, window_type='blackman')
            
        # Let's start with the easy ones!
        self.max = self.getMax()
        
        self.min = self.getMin()
        
        self.avg_height = np.array( self.signals.values() ).mean()
        
        # Let's go with the function fitting
        xdata = np.array( [x for x in sorted(self.signals.keys())] )
        ydata = np.array( [self.signals[x] for x in xdata] )
        plateau, slope, inflection, y0 = fitData(xdata, ydata)
        
        from scipy.integrate import quad
        from DuctApe.Phenome import fitting
        # Some esoteric method for area calculation (to be checked!)
        if plateau:
            area = quad(fitting.gompertz, xdata.min(), xdata.max(), args=(plateau, slope, inflection, y0))
        else:
            area = 0

        print '\t'.join( [self.plate_id, self.well_id] + [str(x) for x in [area, plateau, slope, inflection, y0]] )

class BiologParser(CommonThread):
    '''
    Class BiologParser
    Takes a csv output file and returns a series of objects
    results [PlateCarrier] --> well_id --> BiologRaw
    There could be more than one organism for each biolog file!
    '''
    _statusDesc = {0:'Not started',
               1:'Parsing'}
    
    _substatuses = []
    
    _start = 'Data File'
    _plate = 'Plate Type'
    _strainType = 'Strain Type'
    _sample = 'Sample Number'
    _strainName = 'Strain Name'
    _strainNumber = 'Strain Number'
    _other = 'Other'
    _dataStart = 'Hour'
    
    def __init__(self,infile,queue=Queue.Queue()):
        CommonThread.__init__(self,queue)
        # Biolog
        self.file = infile
        # Results
        self.results = []
        
    def parse(self):
        plate = None
        data = False
        wells = []
        
        tblreader = csv.reader(open(self.file, 'rbU'), delimiter=',')
        for line in tblreader:
            if len(line) == 0:
                continue
            elif self._start in line[0].strip():
                # Do we have to save the old plate?
                if plate:
                    self.results.append(plate)
                data = False
                wells = []
                plate = PlateCarrier()
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
            elif self._other in line[0].strip():
                plate.other = line[1].strip()
            elif self._dataStart in line[0].strip():
                data = True
                for i in range(len(line)):
                    if i == 0:continue
                    x = line[i]
                    if x == '':continue
                    plate.data[x.strip()] = BiologRaw(plate.plate_id, x.strip())
                    plate._idx[i] = x.strip()
                    wells.append(x.strip())
            elif data:
                time = float(line[0])
                for i in range(len(line)):
                    if i == 0:continue
                    x = line[i]
                    if x == '':continue
                    well = plate._idx[i]
                    plate.data[well].addSignal(time, float(x))
        
        return True
    
    def run(self):
        self.updateStatus()
        
        if not self.parse():
            self.sendFailure('Could not parse Biolog file!')
            return

class BiologZero(CommonThread):
    '''
    Class BiologZero
    Takes a list of PlateCarrier objects and subtract the zero signal
    Two modes: 
        - normal (subtract the first well)
        - blank plate (subtract the zero plate --> zero list)
    results [PlateCarrier] --> well_id --> BiologRaw
    '''
    _statusDesc = {0:'Not started',
               1:'Zero subtraction'}
    
    _substatuses = [1]
    
    def __init__(self, data, blank=False, blankData=[], 
                 forceZero = True,
                 queue=Queue.Queue()):
        CommonThread.__init__(self,queue)
        # Biolog
        self.data = data
        self.blank = bool(blank)
        self.blankData = blankData
        self.forceZero = bool(forceZero)
        # Results
        self.results = []
    
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
                        logging.debug('Time %f present in blank plate was not'%(hour)+ 
                                        ' found on plate %s, signal was forced to'%(plate.plate_id)+
                                        ' zero')
                # Reset those hours that are not present in the blank plate
                for hour in sorted(plate.data[well].signals.keys()):
                    if hour not in zplate.data[well].signals:
                        logging.debug('Time %f present in plate %s was not'%(hour, plate.plate_id)+
                                        ' found on blank plate, signal was forced to'+
                                        ' zero')
                        plate.data[well].signals[hour] = 0.1
            
        if not found:
            logging.warning('Blank plate zero subtraction: could not find'+
                                'blank plate %s'%plate.plate_id)
    
    def zeroSubTract(self):
        '''
        Zero subtraction
        '''
        self._maxsubstatus = len(self.data)
        
        if self.blank:
            logging.info('Zero subtraction with blank plate')
        else:
            logging.info('Normal zero subtraction')
        
        for plate in self.data:
            self._substatus += 1
            self.updateStatus(sub=True)
            
            if self.blank:
                self._zeroBlank(plate)
            else:
                self._zeroNormal(plate)
                
        self.results = self.data
                
        return True
    
    def run(self):
        self.updateStatus()
        if not self.zeroSubTract():
            self.sendFailure('Zero subtraction failure!')
            return
        self.resetSubStatus()
        
class BiologPlot(CommonThread):
    '''
    Class BiologPlot
    Takes a list of PlateCarrier objects and creates some plots
    Can work in parallelization
    '''
    _statusDesc = {0:'Not started',
                1:'Making room',
                2:'Preparing data',
                3:'Preparing plots',
                4:'Creating plates plots'}
    
    _substatuses = [2,4]
    
    def __init__(self, data,
                 expname = 'exp', 
                 plateNames = {}, wellNames = {},
                 colors = {}, smooth = True, window = 11,
                 compress = 0,
                 maxsig = None, queue=Queue.Queue()):
        CommonThread.__init__(self,queue)
        # Biolog
        self.data = data
        
        self.expname = expname
        
        # Plot parameters
        self.colors = colors
        self.plateNames = plateNames
        self.wellNames = wellNames
        self.smooth = bool(smooth)
        self.window = int(window)
        self.compress = int(compress)
        self.maxsig = float(maxsig)
        
        # Results
        # PLate_id --> PlotCarrier
        self.results = {}
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
            os.mkdir(path)
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
            logging.warning('Plate %s was not found!'%plate_id)
            return False
        if well_id not in self.results[plate_id].wells:
            logging.warning('Well %s was not found!'%well_id)
            return False
        
        if self.well:
            self.well.clf()
        
        self.well = self.results[plate_id].plotWell(well_id)
        
        return True
    
    def run(self):
        self.updateStatus()
        self.makeRoom()
        
        self._maxsubstatus = len(self.data)
        self.updateStatus()
        for plate in self.data:
            self._substatus += 1
            self.updateStatus(True)
            if plate.plate_id in self.plateNames:
                plate_name = self.plateNames[plate.plate_id]
            else:
                plate_name = ''
            if plate.plate_id not in self.results:
                self.results[plate.plate_id] = PlotCarrier(plate.plate_id, 
                                                    plate_name = plate_name,
                                                    smooth = self.smooth,
                                                    window = self.window,
                                                    compress = self.compress,
                                                    maxsig = self.maxsig)
            
            self.results[plate.plate_id].addData(plate.strain, plate)
            # Sanity check
            if plate.strain not in self.colors:
                logging.error('Color code for strain %s is missing'%plate.strain)
                return
            #
            self.results[plate.plate_id].setColor(plate.strain,
                                                  self.colors[plate.strain])
        self.resetSubStatus()
        
        self.updateStatus()
        for plate_id in self.results:
            self.results[plate_id].preparePlot()
        self.resetSubStatus()
        
        # TODO: more precise here
        # TODO: speed up by implementing threading
        self._maxsubstatus = len(self.results)*96
        self.updateStatus()
        for plate_id in self.results:
            for i in self.results[plate_id].plotAll():
                self._substatus += 1
                self.updateStatus(True)
            # TODO: remember qgraphicspixmapitem for GUI clickable!
            fname = os.path.join(self._room,'%s.png'%plate_id)
            self.results[plate_id].figure.savefig(fname, dpi=150)
            self.results[plate_id].figure.clf()
        self.resetSubStatus()
        
        