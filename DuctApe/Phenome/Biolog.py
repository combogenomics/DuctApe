#!/usr/bin/env python
"""
Biolog

Pheenome library

Classes to handle Biolog data
"""
from DuctApe.Common.CommonThread import CommonThread
from DuctApe.Common.utils import smooth
import Queue
import copy
import csv
import logging
import matplotlib.pyplot as plt

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

class PlotCarrier(object):
    '''
    Class PlotCarrier
    Contains all the distinct strains data for a particular plate
    There can be more than one replica for each strain
    An averaged plate can be added
    '''
    def __init__(self, plate_id, smooth = True, window = 11):
        self.plate_id = plate_id
        self.strains = {}
        self.colors = {}
        
        self.times = None
        
        self.smooth = bool(smooth)
        self.window = int(window)
    
    def _plot(self, dWell):
        '''
        Smooths the signal and plots it
        If there are more than one exp for a strain, the intersection is plotted
        '''
        # Figure creation
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for strain,signals in dWell.iteritems():
            if len(signals) > 1:
                # Intersect!
                logging.warning('Nothing to see here, still not implemented!')
            else:
                # Single plot!
                if self.smooth:
                    signal = smooth(signals[0], window_len=self.window)
                else:
                    signal = signals[0]
                ax.plot(self.times, signal, color=self.colors[strain])
        
        # TODO: Use an unique figure with many subplots?
        return fig
    
    def plot(self):
        '''
        Plot a series of strains
        Returns a series of matplotlib figures (or plots)
        '''
        figures = []
        
        # Check colors
        for strain in self.strains:
            if strain not in self.colors:
                logging.error('Color code for strain %s is missing!'%strain)
                return figures
                
        # Check time concordance
        times = []
        for strain, plates in self.strains.iteritems():
            for plate in plates:
                for well_id, data in plate.data.iteritems():
                    for time in data.signals.keys():
                        if time not in times:
                            times.append(time)
        times.sort()
        self.times = times
        
        # Get also each plate/well pair
        wells = []
        for strain, plates in self.strains.iteritems():
            for plate in plates:
                for well_id, data in plate.data.iteritems():
                    data.fillMissing(self.times)
                    if well_id not in wells:
                        wells.append(well_id) 
        
        # Cycle over each well
        for well_id in wells:
            strain_signals = {}
            for strain, plates in self.strains.iteritems():
                strain_signals[strain] = []
                try:
                    for plate in plates:
                        strain_signals[strain].append( 
                                [plate.data[well_id].signals[hour] for hour in self.times]
                                )
                except:
                    # TODO: add a warning here
                    pass
            
            # Smooth & Plot
            figure = self._plot(strain_signals)
            # TODO: Use a dictionary instead of a list 
            figures.append(figure)
        
        return figures
        
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
        
    def addSignal(self,time,signal):
        self.signals[time] = signal
        
    def fillMissing(self, times):
        '''
        Given a times list, fills the missing values with None
        '''
        for hour in sorted(times):
            if hour not in self.signals:
                self.signals[hour] = None

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
    
    def __init__(self,data, blank=False, blankData=[], queue=Queue.Queue()):
        CommonThread.__init__(self,queue)
        # Biolog
        self.data = copy.deepcopy(data)
        self.blank = bool(blank)
        self.blankData = blankData
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
                        plate.data[well].signals[hour] = 0
            
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