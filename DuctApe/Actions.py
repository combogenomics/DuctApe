#!/usr/bin/env python
"""
Actions

DuctApe Library

All the actions required for the analysis
"""
from DuctApe.Common.Database import DBBase, Project, Genome, Organism
import logging
import os

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('Actions')

################################################################################
# Methods

def dInit(project, wdir='.', name='', descr=''):
    '''
    Initializes a project
    '''
    if not os.path.exists(project):
        create = DBBase(project)
        create.create()
        proj = Project(project)
        tmp = os.path.join(wdir, 'tmp')
        proj.addProject(name=name, description=descr, tmp=tmp)
        logger.info('Project successfully created')
        return True
    else:
        projshort = os.path.split(project)[1]
        logger.warning('Project %s is already present in %s'%(projshort,wdir))
        return False
    
def dRemove(project):
    '''
    Completely removes a project
    '''
    try:
        os.remove(project)
        return True
    except:
        logger.error('Could not remove project %s'%project)
        return False

def dGenomeAdd(project, orgID, filename, name='', descr=''):
    '''
    Add a single genome
    '''
    if not os.path.exists(filename):
        logger.error('Fasta file %s may not be present'%(filename))
        return False
    else:
        org = Organism(project)
        org.addOrg(orgID, name=name, description=descr, orgfile=filename)
        gen = Genome(project)
        gen.addProteome(orgID, filename)
        return True

def dGenomeRemove(project, organisms):
    '''
    Remove all the genomic data about specific organism ID(s)
    '''
    gen = Genome(project)
    for org in organisms:
        gen.delProteome(org)
    return True

def dGenomeClear(project):
    '''
    Clear the genomic tables
    '''
    gen = Genome(project)
    gen.clearAllGenome()
    return True

def dGenomeDirAdd(project, folder):
    if not os.path.exists(folder):
        logger.error('Fasta folder %s may not be present'%(folder))
        return False
    else:
        added = 0
        for infile in os.listdir(folder):
            orgID = infile.split('.')[0]
            filename = os.path.join(folder, infile)
            if os.path.isdir(filename):
                continue
            if not dGenomeAdd(project, orgID, filename):
                logger.error('Could not add genome %s'%infile)
                return False
            added += 1
        if added > 0:
            logger.info('Added %d genomes from %s'%
                    (added, folder))
        else:
            logger.warning('No genomes were added from %s'%folder)
        return True

def isProject(project):
    '''
    Checks if the project file is there
    '''
    if not os.path.exists(project):
        projshort = os.path.split(project)[1]
        logger.warning('Project %s has not been initialized yet!'%projshort)
        return False
    else:
        return True
    
def touchProject(project):
    '''
    Check and update the project file
    '''
    if not isProject(project):
        return False
    else:
        proj = Project(project)
        proj.updateLast()
        logger.debug('%s'%str(proj))
        return True
