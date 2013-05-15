#!/usr/bin/env python
"""
Map2KO

Genome library

Handle a KO search on a local machine (Blast-BBH) or online (KAAS)
"""
from Bio import SeqIO
from ductape.common.commonmultiprocess import CommonMultiProcess
from ductape.common.utils import slice_it
from ductape.genome.blast import Blaster, RunBBH
import Queue
import logging
import os
import webbrowser

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('ductape.map2KO')

################################################################################
# Classes

class KOBBH(object):
    pass

class LocalSearch(CommonMultiProcess):
    '''
    Class localSearch
    '''
    _statusDesc = {0:'Not started',
               1:'Making room', 
               2:'Creating Blast DB',
               3:'Running Blast',
               4:'Running Blast on short proteins',
               5:'Parsing Blast',
               6:'Running BBH',
               7:'Cleaning up'}
    
    _substatuses = [3,4,6]
    
    def __init__(self,query,target,
                 ncpus=1,evalue=1e-50,
                 buildDB=True,bbh=True,recover=False,queue=Queue.Queue()):
        CommonMultiProcess.__init__(self,ncpus,queue)
        # Blast
        self.query = query
        if buildDB:
            self.target = target
            self.db = None
        else:
            self.target = None
            self.db = target
        self.out = []
        self.evalue = float(evalue)
        self.bbh = bool(bbh)
        self.recover = recover
        self.ncpus = int(ncpus)
        self._kohits = []
        self.results = {}
        self._keggroom = None
        self._blast = Blaster()
        
    def makeRoom(self,location=''):
        '''
        Creates a tmp directory in the desired location
        '''
        try:
            path = os.path.abspath(location)
            path = os.path.join(path, 'tmp')
            try:os.mkdir(path)
            except:pass
            path = os.path.join(path, 'blast')
            self._room = path
            os.mkdir(path)
        except:
            logger.debug('Temporary directory creation failed! %s'
                          %path)
        
        # KEGG database path
        try:
            path = os.path.abspath(location)
            path = os.path.join(path, 'tmp')
            try:os.mkdir(path)
            except:pass
            path = os.path.join(path, 'keggdb')
            self._keggroom = path
            os.mkdir(path)
        except:
            logger.debug('Temporary directory creation failed! %s'
                          %path)
        
    def createDB(self):
        self.db = os.path.join(self._keggroom,'KEGGdb') 
        return self._blast.createDB(self.target, 'prot', self.db)
    
    def runBlast(self,short=False):
        lS = []
        for s in SeqIO.parse(open(self.query),'fasta'):
            if short and len(s) <= 30:
                lS.append(s)
            elif not short:lS.append(s)
        self._maxsubstatus = len(lS)
        for seqs in slice_it(lS,10):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return False
            
            if len(seqs) == 0:
                continue
            self._substatus += len(seqs)
            self.updateStatus(sub=True)
            if short:
                query = os.path.join(self._room,
                         'KEGGshort_%d.faa'%self._substatus)
                out = os.path.join(self._room,
                           'KEGGshort_%d.xml'%self._substatus)
            else:
                query = os.path.join(self._room,
                         'KEGG_%d.faa'%self._substatus)
                out = os.path.join(self._room,'KEGG_%d.xml'%self._substatus)
            self.out.append(out)
            # If recovery, skip the unnecessary scans
            if ( self.recover and os.path.exists(query) and
                 os.path.exists(out)):
                # Last test: can it be parsed?
                try:
                    self._blast.parseBlast(out)
                    for hits in self._blast.getHits(self.evalue):
                        pass
                    logger.debug('Skipping slice %s because has already been done'
                                %query)
                    continue
                except:
                    pass
            oseqs = SeqIO.write(seqs,open(query,'w'),'fasta')
            if oseqs != len(seqs):
                logger.warning('Query splitting error! Expected %d, '+
                                'Printed %d'%(len(seqs),oseqs))
            if short:
                res = self._blast.runBlast(query, self.db, out,
                                 evalue = self.evalue,
                                 ncpus = self.ncpus, task='blastp-short')
            else:
                res = self._blast.runBlast(query, self.db, out,
                                 evalue = self.evalue,
                                 ncpus = self.ncpus)
            if not res:
                return False
        return True
    
    def parseBlast(self):
        for out in self.out:
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return False
        
            self._blast.parseBlast(out)
            # Catch the exceptions if the XML is dirty
            try:
                for hits in self._blast.getHits(self.evalue):
                    for hit in hits:
                        if hit.getKO():
                            self._kohits.append(hit)
            except:
                logger.error('Blast results corrupted for file %s'%out)
                return False
        return True
    
    def runBBH(self):
        # Create a DB of the source genome
        sourceDB = os.path.join(self._room,'SOURCEdb') 
        if not self._blast.createDB(self.query, 'prot', sourceDB):
            logger.error('Could not create source DB %s'%sourceDB)
            return False
        self._maxsubstatus = len(self._kohits)
        
        self.initiateParallel()
        
        for hit in self._kohits:
            uniqueid = self.getUniqueID()
            
            if hit.query_len > 30:
                short = False
            else:
                short = True
            
            # Multi process
            obj = RunBBH('Map2KO',hit.query_id,sourceDB,
                    self.db,None,
                    self.evalue,'BLOSUM62',short,uniqueid,
                    kegg = True, ko_entry = hit.hit, ko_id = hit.getKO())
            self._paralleltasks.put(obj)
            
        # Poison pill to stop the workers
        self.addPoison()
        
        while True:
            while not self._parallelresults.empty():
                if self.killed:
                    logger.debug('Exiting for a kill signal')
                    return False
                
                self._substatus += 1
                self.updateStatus(sub=True)
                
                result = self._parallelresults.get()
                
                if not result[2]:
                    logger.error('An error occurred for BBH!')
                    return False
                if result[1] not in self.results and result[0]:
                    self.results[result[1]] = []
                if result[0] and result[0] not in self.results[result[1]]:
                    self.results[result[1]].append(result[0])
                    
            if self.isTerminated():
                break
            
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return False
            
            self.sleeper.sleep(0.1)
            
        # Get the last messages
        while not self._parallelresults.empty():
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return False
            
            self._substatus += 1
            self.updateStatus(sub=True)
            
            result = self._parallelresults.get()
            
            if not result[2]:
                logger.error('An error occurred for BBH!')
                return False
            if result[1] not in self.results and result[0]:
                self.results[result[1]] = []
            if result[0] and result[0] not in self.results[result[1]]:
                self.results[result[1]].append(result[0])
            
        self.killParallel()            

        return True
            
    def run(self):
        self.updateStatus()
        self.makeRoom()

        if self.killed:
            return

        if not self.db:
            self.updateStatus()
            if not self.createDB():
                self.sendFailure('CreateDB failure')
                return
        else:
            self.updateStatus(send=False)
            
        if self.killed:
            return

        self.updateStatus()
        if not self.runBlast():
            self.sendFailure('RunBlast failure')
            return
        self.resetSubStatus()
        
        if self.killed:
            return
        
        self.updateStatus()
        if not self.runBlast(True):
            self.sendFailure('RunBlast (short) failure')
            return
        self.resetSubStatus()
        
        if self.killed:
            return
        
        self.updateStatus()
        if not self.parseBlast():
            self.sendFailure('ParseBlast failure')
            return
        if len(self._kohits) == 0:
            logger.warning('No protein in %s with homology to KO!'%self.query)
            self.sendFailure('No protein in %s with homology to KO!'%self.query)
            self.cleanUp()
            return
        
        if self.killed:
            return
        
        if self.bbh:
            self.updateStatus()
            if not self.runBBH():
                self.sendFailure('BBH failure')
        else:
            for hit in self._kohits:
                if hit.query_id not in self.results:
                    self.results[hit.query_id] = []
                ko = hit.getKO()
                if ko not in self.results[hit.query_id]:
                    self.results[hit.query_id].append(ko)
            self.updateStatus(send=False)   
        self.resetSubStatus()
        
        if self.killed:
            return
        
        try:
            del self.results[None]
        except:pass
        
        # Only ONE KO for each protein
        for k in self.results:
            self.results[k] = self.results[k][0]
        if len(self.results) == 0:
            logger.warning('No protein in %s with BBH to KO!'%self.query)
            self.sendFailure('No protein in %s with BBH to KO!'%self.query)
            self.cleanUp()
            return
        
        if self.killed:
            return
        
        self.updateStatus()
        self.cleanUp()

class OnlineSearch(object):
    '''
    Online KO search using KAAS
    '''
    def __init__(self,bbh=True):
        self.bbh = bbh
        self.results = {}
        self.url = 'http://www.genome.ad.jp/tools/kaas/'
        
    def getExplanation(self):
        msg = ' '.join(['The only way to access to KAAS annotation service is',
                'your browser.\nGo to '+self.url+
                ', where you can add your fasta file(s).',
                '\nThe analysis will take less than one hour.'])
        if self.bbh:
            msg = '\n'.join([msg, 'Select the BBH option and use the appropriate'+
                            ' organisms list.'])
        else:
            msg = '\n'.join([msg, 'Select the SBH option and use the appropriate'+
                            ' organisms list.'])
        msg = '\n'.join([msg, 'When the analysis is finished save the KO list file.'])
        return msg
    
    def openBrowser(self):
        webbrowser.open(self.url, new=2)
        
    def parseKAAS(self,filein):
        for l in open(filein):
            if l.startswith('#'):continue
            s=l.strip().split('\t')
            if len(s) == 1:continue
            self.results[s[0]] = s[1]