#!/usr/bin/env python
"""
Actions

DuctApe Library

All the actions required for the analysis
"""
# TODO: this part must be handled somewhere else
import matplotlib
matplotlib.use('Agg')
#
from ductape.common.utils import slice_it, rgb_to_hex
from ductape.storage.SQLite.database import DBBase, Project, Genome, Organism, \
    Kegg
import logging
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

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
        filename = os.path.abspath(filename)
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
    oCheck = Organism(project)
    for org in organisms:
        if not oCheck.isOrg(org):
            logger.warning('Genome %s is not present: skipping'%org)
            continue
        gen.delProteome(org)
        logger.info('Successfully removed genome %s'%org)
    return True

def dGenomeClear(project):
    '''
    Clear the genomic tables
    '''
    gen = Genome(project)
    gen.clearAllGenome()
    logger.info('Successfully removed all genomic data')
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

def dPanGenomeAdd(project, orthfile):
    '''
    Add an external pangenome
    '''
    if not os.path.exists(orthfile):
        logger.error('Pangenome file %s may not be present'%(orthfile))
        return False
    else:
        orth = {}
        for l in open(orthfile):
            s = l.strip().split('\t')
            if s[0] not in orth:
                orth[s[0]] = []
            orth[s[0]].append(s[1])
        gen = Genome(project)
        gen.addPanGenome(orth)
        
        logger.info('PanGenome size: %d groups'%len(gen.getPanGenome()))
        logger.info('Core size: %d groups'%gen.getLenCore())
        logger.info('Accessory size: %d groups'%gen.getLenAcc())
        logger.info('Unique size: %d groups'%gen.getLenUni())
        
        return True

def dGenomeSetKind(project):
    '''
    Set the kind of genomic project and return its value
    '''
    proj = Project(project)
    proj.getProject()
    org = Organism(project)
    if org.howManyMutants() > 0:
        logger.info('%d mutants are present'%org.howManyMutants())
        proj.setKind('mutants')
        return 'mutants'
    elif org.howMany() == 1:
        logger.info('Just one organism is present')
        proj.setKind('single')
        return 'single'
    else:
        logger.info('%d organisms are present'%org.howMany())
        proj.setKind('pangenome')
        return 'pangenome'

def dGetGenomeSteps(project):
    '''
    Get the analysis that these genomes deserve
    '''
    proj = Project(project)
    proj.getProject()
    status = proj.genome
    pangenome = bool(proj.pangenome)
    kind = dGenomeSetKind(project)
    if kind == 'mutants':
        if status == 'map2ko':
            return ['map2kegg']
        elif status == 'map2kegg':
            return []
        else:
            return ['map2ko', 'map2kegg']
    elif kind == 'single':
        if status == 'map2ko':
            return ['map2kegg']
        elif status == 'map2kegg':
            return []
        else:
            return ['map2ko', 'map2kegg']
    else:
        steps = []
        if not pangenome:
            steps.append('pangenome')
        if status == 'map2ko':
            steps.append('map2kegg')
        elif status == 'map2kegg':
            pass
        else:
            steps.append('map2ko')
            steps.append('pangenome')
        
        return steps

def getPathsReacts(project):
    kegg = Kegg(project)
    # Get the pathway - reaction links
    paths = {}
    for pR in kegg.getPathReacts():
        if pR.path_id in ['path:rn01100','path:rn01110','path:rn01120']:
            continue
        if pR.path_id not in paths:
            paths[pR.path_id] = []
        paths[pR.path_id].append(pR.re_id)
        
    return paths

def prepareColors(dReacts, colorrange):
    if len(dReacts) == 0:
        return {}
    
    maximum = max([dReacts[x] for x in dReacts])
    hexs = {}
    prev = '#FFFFFF'
    i = 1
    for color in slice_it(colorrange, cols=maximum):
        if len(color) == 0:
            hexs[i] = prev
        else:
            hexs[i] = rgb_to_hex(tuple([int(round(x*255))
                              for x in color[-1][:3]])).upper()
        prev = hexs[i]
        i += 1
    
    return hexs

def createLegend(kind, compounds=False):
    '''
    Create a color scheme legend
    '''
    # TODO: a more centralized color scheme
    fig = plt.figure()
    fname = 'legend.png'
    matrix = np.outer(np.arange(0.33,1,0.01),np.ones(7))
    if kind == 'pangenome':
        if compounds:
            pass
        else:
            ax = fig.add_subplot(131)
            ax.imshow(matrix, cmap=cm.Blues, vmin=0, vmax=1)
            ax.axes.get_xaxis().set_visible(False)
            ax.axes.get_yaxis().set_visible(False)
            ax.axes.get_yaxis().set_ticks([])
            ax.axes.get_xaxis().set_ticks([])
            ax.set_title('Core')
            
            ax = fig.add_subplot(132)
            ax.imshow(matrix, cmap=cm.Greens, vmin=0, vmax=1)
            ax.axes.get_xaxis().set_visible(False)
            ax.axes.get_yaxis().set_visible(False)
            ax.axes.get_yaxis().set_ticks([])
            ax.axes.get_xaxis().set_ticks([])
            ax.set_title('Core and Dispensable')
            
            ax = fig.add_subplot(133)
            ax.imshow(matrix, cmap=cm.Oranges, vmin=0, vmax=1)
            ax.axes.get_xaxis().set_visible(False)
            ax.axes.get_yaxis().set_visible(False)
            ax.axes.get_yaxis().set_ticks([])
            ax.axes.get_xaxis().set_ticks([])
            ax.set_title('Dispensable')
            
            fig.savefig(fname)
            
    elif kind == 'single':
        ax = fig.add_subplot(111)
        ax.imshow(matrix, cmap=cm.Greens, vmin=0, vmax=1)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Reactions')
        
        fig.savefig(fname)

    elif kind == 'mutants':
        ax = fig.add_subplot(131)
        ax.imshow(matrix, cmap=cm.Greens, vmin=0, vmax=1)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Wild-type')
        
        ax = fig.add_subplot(132)
        ax.imshow(matrix, cmap=cm.copper_r, vmin=0, vmax=1)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Wild-type and Mutated')
        
        ax = fig.add_subplot(133)
        ax.imshow(matrix, cmap=cm.Reds, vmin=0, vmax=1)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Mutated')
        
        fig.savefig(fname)
        
    return fname

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