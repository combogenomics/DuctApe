#!/usr/bin/env python
"""
CommonThread

Common library

Thread base structure
"""
import Queue
import logging
import os
import shutil
import threading
import time

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('ductape.commonthread')

################################################################################
# Classes

class Status(object):
    '''
    Class Status
    Gives informations about the run status of a specific thread
    '''
    def __init__(self,status=None,msg=None,maxstatus=None,
                    substatus=None,submsg=None,maxsubstatus=None,
                    fail=False):
        self.status = status
        self.msg = msg
        self.maxstatus = maxstatus
        #
        self.substatus = substatus
        self.submsg = submsg
        self.maxsubstatus = maxsubstatus
        # Fail msg?
        self.fail = fail

class CommonThread(threading.Thread):
    '''
    Class CommonThread: Common operations for a threading class
    '''
    _statusDesc = {0:'Not started',
               1:'Making room', 
               3:'Cleaning up'}
    
    _substatuses = []
    
    def __init__(self,queue=Queue.Queue()):
        threading.Thread.__init__(self)
        # Thread
        self.msg = queue
        self._status = 0
        self._maxstatus = len(self._statusDesc)
        self._substatus = 0
        self._maxsubstatus = 0
        self._room = None
        self.killed = False
        
    def getStatus(self):
        return self._statusDesc[self._status]
    
    def getMaxStatus(self):
        return self._maxstatus
    
    def getMaxSubStatus(self):
        return self._maxsubstatus
    
    def getSubStatuses(self):
        return self._substatuses
    
    def resetSubStatus(self):
        self._substatus = 0
        self._maxsubstatus = 0
        
    def makeRoom(self,location=''):
        '''
        Creates a tmp directory in the desired location
        '''
        try:
            path = os.path.abspath(location)
            path = os.path.join(path, 'tmp')
            self._room = path
            os.mkdir(path)
        except:
            logger.debug('Temporary directory creation failed! %s'
                          %path)
    
    def startCleanUp(self):
        '''
        Removes the temporary directory
        '''
        if os.path.exists(self._room):
            logger.debug('Removing the old results directory (%s)'%
                         self._room)
            shutil.rmtree(self._room, True)
    
    def cleanUp(self):
        '''
        Removes the temporary directory
        '''
        shutil.rmtree(self._room, True)
        
    def run(self):
        self.updateStatus()
        self.makeRoom()

        self.updateStatus()
        self.cleanUp()
            
    def sendFailure(self,detail='Error!'):
        msg = Status(fail=True,
                     msg=detail)
        self.msg.put(msg)
        # Give some time for the message to arrive
        time.sleep(0.1)
        
    def updateStatus(self,sub=False,send=True):
        if not sub:
            self._status += 1
        if not send:
            return
        if self._status in self._substatuses:
            msg = Status(status=self._status,msg=self.getStatus(),
                         maxstatus=self.getMaxStatus(),
                         substatus=self._substatus,
                         maxsubstatus=self.getMaxSubStatus())
        else:
            msg = Status(status=self._status,msg=self.getStatus(),
                         maxstatus=self.getMaxStatus())
        self.msg.put(msg)
        
    def kill(self):
        self.killed = True
