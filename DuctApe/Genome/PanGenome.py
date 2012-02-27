#!/usr/bin/env python
"""
PanGenome

Genome library

Uses a serial-BBH approach to compute a pangenome of the desired organisms list
"""
from Bio import SeqIO
from DuctApe.Common.CommonMultiProcess import CommonMultiProcess
from DuctApe.Genome.Blast import Blaster, RunBBH
import Queue
import logging
import os
import shutil

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('PanGenome')

################################################################################
# Classes

class PanGenomer(CommonMultiProcess):
    '''
    Class panGenomer
    '''
    _statusDesc = {0:'Not started',
               1:'Making room',
               2:'Creating Blast DBs',
               3:'Running Blast BBHs',
               4:'Crafting the PanGenome',
               5:'Cleaning up'}
    
    _substatuses = [2,3]
    
    def __init__(self,organisms,
                 ncpus=1,evalue=1e-50,
                 recover=False,prefix='',
                 matrix='BLOSUM80',queue=Queue.Queue()):
        CommonMultiProcess.__init__(self,ncpus,queue)
        # Blast
        self.organisms = list(organisms)
        self.dbs = {}
        self._prot2orgs = {}
        self.out = []
        self.evalue = float(evalue)
        # TODO: implement recovery
        self.recover = recover
        #
        self.results = {}
        self._blast = Blaster()
        self._pangenomeroom = None
        self.prefix = prefix.rstrip('_')
        self.matrix = matrix
        self._already = []
        # Results
        self.orthologs = {}
        self.core = []
        self.accessory = []
        self.unique = []

    def makeRoom(self,location=''):
        '''
        Creates a tmp directory in the desired location
        '''
        try:
            path = os.path.abspath(location)
            path = os.path.join(path, 'tmp')
            try:os.mkdir(path)
            except:pass
            path = os.path.join(path, 'pangenomeDBs')
            self._room = path
            os.mkdir(path)
        except:
            logger.debug('Temporary directory creation failed! %s'
                          %path)
        
        try:
            path = os.path.abspath(location)
            path = os.path.join(path, 'tmp')
            try:os.mkdir(path)
            except:pass
            path = os.path.join(path, 'pangenome')
            self._pangenomeroom = path
            os.mkdir(path)
        except:
            logger.debug('Temporary directory creation failed! %s'
                          %path)
    
    def cleanUp(self):
        '''
        Removes the temporary directory
        '''
        shutil.rmtree(self._room, True)
        shutil.rmtree(self._pangenomeroom, True)
                
    def createDB(self):
        dbindex = 0
        
        self._maxsubstatus = len(self.organisms)
        
        for org in self.organisms:
            self._substatus += 1
            self.updateStatus(sub=True)
            
            seqs = [seq.id for seq in SeqIO.parse(open(org),'fasta')]
            for seqid in seqs:
                if seqid in self._prot2orgs:
                    logging.warning('Protein %s present as duplicate!'%seqid)
                    return False
                self._prot2orgs[seqid] = org            
            self.dbs[org] = os.path.join(self._room,str(dbindex)) 
            res = self._blast.createDB(org, 'prot', self.dbs[org])
            if not res:
                logger.error('Could not create DB for %s'%org)
                return False
            dbindex += 1
        return True
    
    def serialBBH(self):
        orthindex = 1
        
        self._maxsubstatus = len(self._prot2orgs)
        
        for org in self.organisms:
            seqs = [seq for seq in SeqIO.parse(open(org),'fasta')]
            # Iterate over each protein
            for seq in seqs:
                self._substatus += 1
                self.updateStatus(sub=True)
                
                if seq.id in self._already:
                    continue
                orthname = self.prefix + str(orthindex)
                orgsincluded = [org]
                self.orthologs[orthname] = [seq.id]
                query = os.path.join(self._pangenomeroom,str(self._substatus))
                if SeqIO.write([seq], open(query,'w'), 'fasta') == 0:
                    logger.error('Error writing sequence %s to file'%seq.id)
                    return False
                
                self.initiateParallel()
                
                # Iterate over each other organism
                for otherorg in self.organisms:
                    if org == otherorg:
                        continue
                    # Go fot it!
                    if len(seq) < 30:
                        short = True
                    else:
                        short = False
                    
                    uniqueid = self.getUniqueID()
                    
                    # Multi process
                    obj = RunBBH(query,seq.id,self.dbs[org],
                            self.dbs[otherorg],otherorg,
                            self.evalue,self.matrix,short,uniqueid)
                    self._paralleltasks.put(obj)
                    
                # Poison pill to stop the workers
                self.addPoison()
                
                while True:
                    while not self._parallelresults.empty():
                        result = self._parallelresults.get()
                        
                        if not result[2]:
                            logger.error('An error occurred for BBH on query %s'%seq.id+
                                         ' and target %s'%result[1])
                            return False
                        if result[0] and result[0] not in self._already:
                            self.orthologs[orthname].append(result[0])
                            orgsincluded.append(result[1])
                            self._already.append(result[0])
                            
                    if self.isTerminated():
                        break
                    
                    self.sleeper.sleep(0.1)
                    
                while not self._parallelresults.empty():
                    result = self._parallelresults.get()
                    
                    if not result[2]:
                        logger.error('An error occurred for BBH on query %s'%seq.id+
                                     ' and target %s'%result[1])
                        return False
                    if result[0] and result[0] not in self._already:
                        self.orthologs[orthname].append(result[0])
                        orgsincluded.append(result[1])
                        self._already.append(result[0])
                
                self.killParallel()
                
                if len(orgsincluded) < len(self.organisms):
                    logging.debug('Additional search on missing organisms for'+
                                  ' ortholog %s'%orthname)
                    for otherprotein in self.orthologs[orthname]:
                        if otherprotein == seq.id:
                            continue
                        neworg = self._prot2orgs[otherprotein]
                        bFound = False
                        for seq in SeqIO.parse(open(neworg),'fasta'):
                            if seq.id == otherprotein:
                                bFound = True
                                break
                        if not bFound:
                            logger.error('%s not found!'%otherprotein)
                            return False
                        if SeqIO.write([seq], open(query,'w'), 'fasta') == 0:
                            logger.error('Error writing sequence %s to file'
                                         %seq.id)
                            return False
                        
                        self.initiateParallel()
                        
                        for evenneworg in self.organisms:
                            if evenneworg in orgsincluded:
                                continue
                            # Go fot it!
                            if len(seq) < 30:
                                short = True
                            else:
                                short = False
                            
                            uniqueid = self.getUniqueID()
                    
                            # Multi process
                            obj = RunBBH(query,otherprotein,self.dbs[neworg],
                                    self.dbs[evenneworg],evenneworg,
                                    self.evalue,self.matrix,short,uniqueid)
                            self._paralleltasks.put(obj)
                            
                        # Poison pill to stop the workers
                        self.addPoison()
                        
                        while True:
                            while not self._parallelresults.empty():
                                result = self._parallelresults.get()
                                
                                if not result[2]:
                                    logger.error('An error occurred for BBH on query %s'%seq.id+
                                                 ' and target %s'%result[1])
                                    return False
                                if result[0] and result[0] not in self._already:
                                    self.orthologs[orthname].append(result[0])
                                    orgsincluded.append(result[1])
                                    self._already.append(result[0])
                            
                            if self.isTerminated():
                                break
                            
                            self.sleeper.sleep(0.1)
                        
                        while not self._parallelresults.empty():
                            result = self._parallelresults.get()
                            
                            if not result[2]:
                                logger.error('An error occurred for BBH on query %s'%seq.id+
                                             ' and target %s'%result[1])
                                return False
                            if result[0] and result[0] not in self._already:
                                self.orthologs[orthname].append(result[0])
                                orgsincluded.append(result[1])
                                self._already.append(result[0])
                        
                        self.killParallel()
                
                os.remove(query)        
                orthindex += 1
        return True
    
    def packPanGenome(self):
        for g in self.orthologs:
            if len(self.orthologs[g]) == len(self.organisms):
                self.core.append(g)
            elif len(self.orthologs[g]) == 1:
                self.unique.append(g)
            else:
                self.accessory.append(g)
    
    def run(self):
        self.updateStatus()
        self.makeRoom()
        
        self.updateStatus()
        if not self.createDB():
            self.sendFailure('Create DBs failure!')
            self.cleanUp()
            return
        self.resetSubStatus()
            
        self.updateStatus()
        if not self.serialBBH():
            self.sendFailure('Serial BBH failure!')
            self.killParallel()
            self.cleanUp()
            return
        self.resetSubStatus()
        
        self.updateStatus()
        self.packPanGenome()
        
        self.updateStatus()
        self.cleanUp()
