#!/usr/bin/env python
"""
Actions Parallel

DuctApe Library

All the actions required for the analysis (parallel utilitites)
"""
from ductape.common.commonmultiprocess import CommonMultiProcess
from ductape.storage.SQLite.database import Kegg
import logging
import Queue

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('ductape.actionsparallel')

################################################################################
# Classes

class RunPath(object):
    def __init__(self, project, path_id):
        self._project = project
        self.path_id = path_id
    
    def __call__(self):
        db = Kegg(self._project)
        allr = [r for r in db.getMappedRPairsReact(self.path_id)]
        ecore, edisp, eacc, euni = db.getExclusiveRPairsReact(self.path_id)
        dpangenome = {'all':allr,
                        'core':ecore, 'dispensable':edisp,
                          'accessory':eacc, 'unique':euni}

        return self.path_id, dpangenome

class PathPanGenomer(CommonMultiProcess):
    '''
    Class PathPangenomer
    '''
    _statusDesc = {0:'Not started',
               1:'Analyzing pangenomic pathways'}
    
    _substatuses = [1]
    
    def __init__(self,project, paths,
                 ncpus=1,queue=Queue.Queue()):
        CommonMultiProcess.__init__(self,ncpus,queue)
        
        # DB name
        self._project = project
        
        # To-be-analyzed pathways
        self.paths = paths
        
        # path_id --> dpangenome
        self.result = {}
    
    def analyzePaths(self):
        self._maxsubstatus = len(self.paths)
        
        self.initiateParallel()
        for path in self.paths:
            # Multi process
            obj = RunPath(self._project, path)
            self._paralleltasks.put(obj)
                
        # Poison pill to stop the workers
        self.addPoison()
        
        while True:
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
                 
            while not self._parallelresults.empty():
                if self.killed:
                    logger.debug('Exiting for a kill signal')
                    return
                
                path_id, dpangenome = self._parallelresults.get()
                self.result[path_id] = dpangenome
                
                self._substatus += 1
                self.updateStatus(sub=True)
                    
            if self.isTerminated():
                break
            
            self.sleeper.sleep(0.1)
            
        while not self._parallelresults.empty():
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            path_id, dpangenome = self._parallelresults.get()
            self.result[path_id] = dpangenome
            
            self._substatus += 1
            self.updateStatus(sub=True)
        
        self.killParallel()
        
        return True
    
    def run(self):
        self.updateStatus()
        if not self.analyzePaths():
            self.sendFailure('Could not analyze pathways!')
            return
        self.resetSubStatus()
