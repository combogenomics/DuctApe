#!/usr/bin/env python
"""
Actions for terminal

DuctApe Library

All the actions required for the analysis (terminal utilitites)
"""
from ductape.storage.SQLite.database import Project, Kegg
import logging

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('ductape.actionsterm')

################################################################################
# Methods

def fetchKegg(project, keeptrying=False):
    from ductape.kegg.kegg import KeggNet, KeggAPI, BaseKegg
    from ductape.terminal import RunThread
    
    # Check if we have to fetch the whole kegg DB
    fetch = False
    proj = Project(project)
    
    logger.info('Checking connectivity')
    bk = BaseKegg()
    try:
        bk.checkConnection()
    except Exception as e:
        logger.error(str(e))
        return False
    
    k = KeggAPI()
    try:
        k.getDBVersion()
        release = k.result[1]
    except Exception as e:
        logger.warning('Could not fetch the KEGG DB version (%s)'%str(e))
        release = None
    if proj.isKegg():
        if release and proj.kegg < release:
            logger.warning('A new KEGG DB version is available (%s, was %s)'%
                           (str(release), str(proj.kegg)))
            fetch = True
    else:
        fetch = True
  
    if fetch:
        logger.info('Fetching the whole KEGG metabolic map')
        if release:
            logger.info('KEGG DB release %s'%str(release))
        kegg = Kegg(project)
        
        knet = KeggNet(keeptrying=keeptrying)
        if not RunThread(knet):
            return False
        
        # Details
        kegg.addPathways(knet.result.path)
        logger.info('Added %d Path IDs'%len(knet.result.path))
        kegg.addReactions(knet.result.react)
        logger.info('Added %d Re IDs'%len(knet.result.react))
        kegg.addCompounds(knet.result.comp)
        logger.info('Added %d Co IDs'%len(knet.result.comp))
        kegg.addRPairs(knet.result.rpair)
        logger.info('Added %d RPair IDs'%len(knet.result.rpair))
        # Links
        kegg.addPathReacts(knet.result.pathreact)
        kegg.addReactComps(knet.result.reactcomp)
        kegg.addCompReacts(knet.result.compreact)
        kegg.addPathComps(knet.result.pathcomp)
        kegg.addReactRPairs(knet.result.reactrpair)
        kegg.addRPairReacts(knet.result.rpairreact)
        logger.info('Added Kegg links')
        # HTML maps
        kegg.addPathHtml(knet.result.pathmaps)
        logger.info('Added Kegg maps')
        
        # Add the release version
        if release:
            proj.setKegg(release)
    else:
        logger.info('KEGG db is up-to-date')
        
    return True
