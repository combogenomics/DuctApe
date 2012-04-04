#!/usr/bin/env python
"""
Terminal

DuctApe Library

Utilities for command line DuctApe programs
"""
from DuctApe.Common.TerminalProgress import ProgressBar, TerminalController
import logging
import sys
import time

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('Terminal')

################################################################################
# Methods

def RunThread(obj):
    obj.start()
    
    prg = None
    substatus = 0
    sub = False
    while True:
        try:
            time.sleep(0.5)
            while not obj.msg.empty():
                msg = obj.msg.get()
                if msg:
                    if msg.status in obj.getSubStatuses():
                        if prg is None or substatus != msg.status:
                            substatus = msg.status
                            sub = True
                            logger.info('%s - %s'%(msg.status, msg.msg))
                            prg = ProgressBar(TerminalController(sys.stdout), msg.msg)
                    if msg.fail:
                        logger.error('Failure! %s'%msg.msg)
                        return False
                    elif not msg.substatus and not sub:
                        logger.info('%s - %s'%(msg.status, msg.msg))
                        prg = None
                    else:
                        if msg.maxsubstatus != 0:
                            prg.update(msg.substatus/float(msg.maxsubstatus),
                                'Item %d on %d total'%(msg.substatus, msg.maxsubstatus))
                    sub = False
                    
            if not obj.isAlive():
                break
        except KeyboardInterrupt:
            logger.warning('Got keyboard interruption -- Aborting')
            obj.kill()
            obj.join()
            logger.warning('Got keyboard interruption -- Aborted')
            return False
        
    return True
