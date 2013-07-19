#!/usr/bin/env python
"""
Terminal

DuctApe Library

Utilities for command line DuctApe programs
"""
from ductape.common.terminalprogress import ProgressBar, TerminalController
import logging
import sys
import time

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('ductape.terminal')

################################################################################
# Methods

def RunThread(obj):
    # Check terminal log state: if lower then INFO
    # Do not display the progressbars
    progress = True
    for hand in logger.parent.handlers:
        if type(hand) == logging.StreamHandler:
            if hand.level < logging.INFO:
                progress = False
                break
    
    obj.start()
    
    prg = None
    substatus = 0
    sub = False
    while True:
        try:
            time.sleep(0.01)
            while not obj.msg.empty():
                msg = obj.msg.get()
                if msg:
                    if msg.status in obj.getSubStatuses():
                        if prg is None or substatus != msg.status:
                            substatus = msg.status
                            sub = True
                            logger.info('%s - %s'%(msg.status, msg.msg))
                            if progress:
                                prg = ProgressBar(TerminalController(sys.stdout), msg.msg)
                            else:
                                prg = ''
                    if msg.fail:
                        logger.error('Failure! %s'%msg.msg)
                        return False
                    elif not msg.substatus and not sub:
                        logger.info('%s - %s'%(msg.status, msg.msg))
                        prg = None
                    else:
                        if msg.maxsubstatus != 0 and progress:
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
