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
        logger.info('Added genome %s, having %d proteins'%
                    (orgID, gen.howMany(orgID)))
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
    '''
    Add a series of genomes contained in a directory
    '''
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

def dGenomeMutAdd(project, mutID, mutparent, mutfasta, kind, name='', descr=''):
    '''
    Check and add a mutant
    '''
    if not os.path.exists(mutfasta):
        logger.error('Fasta file %s may not be present'%(mutfasta))
        return False
    else:
        org = Organism(project)
        if not org.isOrg(mutparent):
            logger.error('Parent organism %s not present!'%mutparent)
            return False
        elif org.isMutant(mutparent):
            logger.error('Parent organism %s cannot be a mutant!'%mutparent)
            return False
        parents = len(org) - org.howManyMutants()
        if parents != 1:
            logger.error('Only one parent is allowed!')
            return False
        org.addOrg(mutID, name=name, description=descr, orgfile=mutfasta,
                   mutant=True, reference=mutparent, mkind=kind)
        gen = Genome(project)
        gen.addProteome(mutID, mutfasta)
        logger.info('Mutant %s (%s) added, having %d mutated genes'
                    %(mutID, org.getOrg(mutID).mkind,gen.howMany(mutID)))
        return True
    
def dGetGenomeSteps(project):
    '''
    Get the analysis that these genomes deserve
    '''
    proj = Project(project)
    proj.getProject()
    status = proj.genome
    org = Organism(project)
    gen = Genome(project)
    if org.howManyMutants() > 0:
        logger.info('%d mutants are present'%org.howManyMutants())
        proj.setKind('mutants')
        if status == 'map2ko':
            return ['map2kegg']
        elif status == 'map2kegg':
            return []
        else:
            return ['map2ko', 'map2kegg']
    elif org.howMany() == 1:
        logger.info('Just one organism is present')
        proj.setKind('single')
        if status == 'map2ko':
            return ['map2kegg']
        elif status == 'map2kegg':
            return []
        else:
            return ['map2ko', 'map2kegg']
    else:
        logger.info('%d organisms are present'%org.howMany())
        proj.setKind('pangenome')
        if status == 'pangenome':
            return ['map2ko', 'map2kegg']
        elif status == 'map2ko':
            if len(gen.getPanGenome()) == 0:
                return ['pangenome', 'map2kegg']
            else:
                return ['map2kegg']
        elif status == 'map2kegg':
            return []
        else:
            return ['pangenome', 'map2ko', 'map2kegg']
    
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

def prepareDir(wdir, tdir):
    '''
    Prepare the temp directory
    '''
    tmp = os.path.join(wdir,tdir)
    if os.path.exists(tmp):
        return True
    else:
        try:
            os.mkdir(wdir)
        except:
            pass
        try:
            os.mkdir(tmp)
        except:
            logger.error('Could not create tmp directory %s'%tmp)
            return False
        return True