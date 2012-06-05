#!/usr/bin/env python
"""
CommonMultiProcess

Common library

MultiProcess base structure
"""
from ductape.common.commonthread import CommonThread
from multiprocessing.queues import Queue
import logging
import multiprocessing
import time

# Consumer borrowed from http://broadcast.oreilly.com/
# EINTR fix borrowed from Boyd Waters

logger = logging.getLogger('ductape.commonmultiprocess')

class SafeSleep(object):
    '''
    IOError safe sleep
    '''
    def sleep(self,seconds):
        '''
        Sleeps for a certain amount of seconds
        Raises an exception if too many errors are encountered
        '''
        dt = 1e-3
        while dt < 1:
            try:
                time.sleep(seconds)
                return
            except IOError:
                logger.warning('IOError encountered in SafeSleep sleep()')
                try:
                    time.sleep(dt)
                except:pass
                dt *= 2
                
        e = IOError('Unrecoverable error')
        raise e

class SafeQueue(Queue):
    '''
    IOError safe multiprocessing Queue
    '''
    def __init__(self):
        Queue.__init__(self)
        
    def empty(self):
        '''
        Returns True if the Queue is empty, False otherwise
        Raises an exception if too many errors are encountered 
        '''
        dt = 1e-3
        while dt < 1:
            try:
                isEmpty = Queue.empty(self)
                return isEmpty
            except IOError:
                logger.warning('IOError encountered in SafeQueue empty()')
                try:
                    time.sleep(dt)
                except:pass
                dt *= 2
                
        e = IOError('Unrecoverable error')
        raise e

    def get(self):
        '''
        Get the element in the queue
        Raises an exception if it's empty or if too many errors are
        encountered
        '''
        dt = 1e-3
        while dt < 1:
            try:
                element = Queue.get(self)
                return element
            except IOError:
                logger.warning('IOError encountered in SafeQueue get()')
                try:
                    time.sleep(dt)
                except:pass
                dt *= 2
                
        e = IOError('Unrecoverable error')
        raise e
    
    def put(self,element):
        '''
        Put the element in the queue
        Raises an exception if too many errors are
        encountered
        '''
        dt = 1e-3
        while dt < 1:
            try:
                Queue.put(self,element)
                return
            except IOError:
                logger.warning('IOError encountered in SafeQueue put()')
                try:
                    time.sleep(dt)
                except:pass
                dt *= 2
                
        e = IOError('Unrecoverable error')
        raise e

class Consumer(multiprocessing.Process):
    
    def __init__(self, 
                 task_queue = multiprocessing.Queue(),
                 result_queue = multiprocessing.Queue()):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.sleeper = SafeSleep()

    def run(self):
        while True:
            next_task = self.task_queue.get()
            self.sleeper.sleep(0.01)
            if next_task is None:
                # Poison pill means we should exit
                break
            answer = next_task()
            self.result_queue.put(answer)
        return
    
class CommonMultiProcess(CommonThread):
    '''
    Class CommonMultiProcess
    A Thread that can perform multiprocesses
    '''
    def __init__(self,ncpus=1, queue=Queue()):
        CommonThread.__init__(self,queue)
        
        self.ncpus = int(ncpus)
        # Parallelization
        self._parallel = None
        self._paralleltasks = SafeQueue()
        self._parallelresults = SafeQueue()
        self.sleeper = SafeSleep()
        
        # ID
        self._unique = 0
        
    def getUniqueID(self):
        self._unique += 1
        return self._unique
    
    def initiateParallel(self):
        self._parallel = [Consumer(self._paralleltasks,self._parallelresults)
                          for x in range(self.ncpus)]
        for consumer in self._parallel:
            consumer.start()
            
    def addPoison(self):
        for consumer in self._parallel:
            self._paralleltasks.put(None)

    def isTerminated(self):
        for consumer in self._parallel:
            if consumer.is_alive():
                return False
        return True

    def killParallel(self):
        for consumer in self._parallel:
            consumer.terminate()
