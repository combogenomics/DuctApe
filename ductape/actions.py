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
from Bio import SeqIO
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

def dGenomeStats(project, doPrint=True):
    # Which project are we talking about?
    kind = dGenomeSetKind(project)
    
    proj = Project(project)
    organism = Organism(project)
    genome = Genome(project)
    kegg = Kegg(project)
    
    if kind == 'single' or kind == 'pangenome':
        logger.info('Single genomes stats')
        # Single genomes stats
        # Header
        header = '\t'.join( ['ID', 'name', 'description', 'proteome size',
                                'mapped to kegg', 'KEGG orthology IDs',
                                'pathways', 'reactions'] )
        if doPrint:
            print header
        else:
            logger.info(header)
        
        lOrg = []
        for org in organism.getAll():
            org_id = org.org_id
            name = org.name if org.name else 'NONE'
            description = org.description if org.description else 'NONE'
            
            prots = genome.howMany(org_id)
            
            mapped, ko, react, path = (kegg.howManyMapped(org_id),
                                        kegg.howManyKO(org_id),
                                        kegg.howManyReactions(org_id),
                                        kegg.howManyPathways(org_id))
            
            stats = '\t'.join( [str(x) for x in [org_id, name, description,
                                                 prots, mapped, ko, path,
                                                 react]] )
            if doPrint:
                print stats
            else:
                logger.info(stats)
                
            lOrg.append([org_id, prots, mapped, react])
            
        plotMapBars(lOrg, 'Single genomes statistics', 'single')
        
        if proj.isPanGenome():
            logger.info('Pangenome stats')
            # Pangenome stats
            # Header
            header = '\t'.join( ['kind', 'size',
                                    'mapped to kegg', 'KEGG orthology IDs',
                                    'pathways', 'reactions'] )
            if doPrint:
                print header
            else:
                logger.info(header)
                
            core, acc, uni = (genome.getLenCore(), genome.getLenAcc(),
                              genome.getLenUni())

            stats = []
            stats.append('\t'.join( [str(x) for x in ['core', core,
                                 kegg.howManyMapped(pangenome='core'),
                                 kegg.howManyKO(pangenome='core'),
                                 kegg.howManyPathways(pangenome='core'),
                                 kegg.howManyReactions(pangenome='core')]]))
            stats.append('\t'.join( [str(x) for x in ['accessory', acc,
                                 kegg.howManyMapped(pangenome='accessory'),
                                 kegg.howManyKO(pangenome='accessory'),
                                 kegg.howManyPathways(pangenome='accessory'),
                                 kegg.howManyReactions(pangenome='accessory')]]))
            stats.append('\t'.join( [str(x) for x in ['unique', uni,
                                 kegg.howManyMapped(pangenome='unique'),
                                 kegg.howManyKO(pangenome='unique'),
                                 kegg.howManyPathways(pangenome='unique'),
                                 kegg.howManyReactions(pangenome='unique')]]))
            
            for stat in stats:
                if doPrint:
                    print stat
                else:
                    logger.info(stat)
            
            lPanGenome = [['Core', core, kegg.howManyMapped(pangenome='core'),
                           kegg.howManyReactions(pangenome='core')],
                          ['Accessory', acc,
                           kegg.howManyMapped(pangenome='accessory'),
                           kegg.howManyReactions(pangenome='accessory')],
                          ['Unique', uni,
                           kegg.howManyMapped(pangenome='unique'),
                           kegg.howManyReactions(pangenome='unique')]]
 
            plotMapBars(lPanGenome, 'PanGenome statistics', 'pangenome_stats')
            plotPanGenome(core, acc, uni)
    
    elif kind == 'mutants':
        refs = [org.org_id
                    for org in organism.getAll()
                    if not organism.isMutant(org.org_id)]
        
        # Header
        header = '\t'.join( ['ID', 'name', 'description', 'kind', 'proteome size',
                                'mapped to kegg', 'reactions'] )
        
        for ref_id in refs:
            logger.info('Mutants of %s stats'%ref_id)
            
            if doPrint:
                print header
            else:
                logger.info(header)
            
            muts = [x for x in organism.getOrgMutants(ref_id)]
            
            lOrg = []
            for org_id in [ref_id] + muts:
                org = organism.getOrg(org_id)
                
                name = org.name if org.name else 'NONE'
                description = org.description if org.description else 'NONE'
                
                mkind = org.mkind if org.mkind in ['deletion', 'insertion'] else 'wild-type'
                
                if mkind not in ['deletion', 'insertion']:
                    prots = genome.howMany(org_id)
                elif mkind == 'deletion':
                    prots = genome.howMany(ref_id) - genome.howMany(org_id)
                elif mkind == 'insertion':
                    prots = genome.howMany(ref_id) + genome.howMany(org_id)
                
                mapped, react = (kegg.howManyMapped(org_id),
                                kegg.howManyReactions(org_id))
        
                if mkind == 'deletion':
                    mapped = kegg.howManyMapped(ref_id) - mapped
                    react = kegg.howManyReactions(ref_id) - react
                elif mkind == 'insertion':
                    mapped += kegg.howManyMapped(ref_id)
                    react += kegg.howManyReactions(ref_id)
                
                stats = '\t'.join( [str(x) for x in [org_id, name, description,
                                                 mkind, prots, mapped,
                                                 react]] )
                if doPrint:
                    print stats
                else:
                    logger.info(stats)
                
                lOrg.append([org_id, prots, mapped, react])
        
            plotMapBars(lOrg, 'Wild-type (%s) and mutants statistics'%ref_id,
                        '%s'%ref_id)
    
    else:
        logger.info('No statistics can be computed at this time')
        return False

    return True
        
def dGenomeExport(project):
    # Is there something to be exported?
    organism = Organism(project)
    
    if organism.howMany() == 0:
        logger.info('No genomic data can be exported at this time')
        return False
    else:
        logger.info('Exporting protein data')
        
        genome = Genome(project)
        
        for org in organism.getAll():
            nprots = SeqIO.write([x for x in genome.getRecords(org.org_id)],
                        open('%s.faa'%org.org_id,'w'), 'fasta')
            logger.info('Saved %d proteins from %s (%s)'%(nprots,
                                                          org.org_id,
                                                          '%s.faa'%org.org_id))
            
        logger.info('Exporting Kegg data')
        
        logger.info('Exporting KO map data')
        
        kegg = Kegg(project)
        
        for org in organism.getAll():
            fname = 'ko_%s.tsv'%org.org_id
            fout = open(fname,'w')
            fout.write('protein_id\tko_id\n')
            i = 0
            for prot_id, ko_id in kegg.getAllKO(org.org_id):
                fout.write('%s\t%s\n'%(prot_id, ko_id))
                i += 1
            fout.close()
            logger.info('Saved %d KO links for %s (%s)'%(i, org.org_id,
                                                         fname))
            
        logger.info('Exporting Kegg reactions data')
        
        for org in organism.getAll():
            fname = 'reactions_%s.tsv'%org.org_id
            fout = open(fname,'w')
            fout.write('protein_id\treaction_id\n')
            i = 0
            for prot_id, re_id in kegg.getAllReactions(org.org_id):
                fout.write('%s\t%s\n'%(prot_id, re_id))
                i += 1
            fout.close()
            logger.info('Saved %d Kegg reactions links for %s (%s)'%
                        (i, org.org_id, fname))
            
        proj = Project(project)
        
        if proj.isPanGenome():
            logger.info('Exporting pangenome data')
            
            fname = 'pangenome.tsv'
            fout = open(fname,'w')
            fout.write('ortholog_id\tprotein_id\n')
            dG = genome.getPanGenome()
            for group, prots in dG.iteritems():
                for prot in prots:
                    fout.write('%s\t%s\n'%(group,prot))
            fout.close()
            
            logger.info('Exported %d orthologs (%s)'%(len(dG),fname))
            
            fname = 'pangenome_category.tsv'
            fout = open(fname,'w')
            fout.write('ortholog_id\tkind\torganisms\n')
            dG = genome.getPanGenomeOrgs()
            for group in genome.getCore():
                fout.write('%s\t%s\t%s\n'%(group.group_id,
                                           'core',
                                           '-'.join(dG[group.group_id])))
            for group in genome.getAcc():
                fout.write('%s\t%s\t%s\n'%(group.group_id,
                                           'accessory',
                                           '-'.join(dG[group.group_id])))
            for group in genome.getUni():
                fout.write('%s\t%s\t%s\n'%(group.group_id,
                                           'unique',
                                           '-'.join(dG[group.group_id])))
            fout.close()
            
            logger.info('Exported orthologs informations (%s)'%fname)
    
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
    elif org.howMany() == 0:
        logger.info('No organisms are present yet')
        return None
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

def plotMapBars(lOrg, title, fname):
    '''
    Plot histograms for Kegg mapping statistics
    '''
    plt.clf()
    space = np.array([0.0, 0.2, 0.4])
    maxprots = max([x[1] for x in lOrg])
    
    for data in lOrg:
        index = float(lOrg.index(data))
        patch = plt.bar(space + index, np.array(data[1:]),
                width=0.2,
                color=['blue', 'orange', 'red'])
        if index == 0:
            patch[0].set_label('Size')
            patch[1].set_label('Mapped to Kegg')
            patch[2].set_label('Kegg reactions')
    
    plt.xticks([0.2 + 1 * x for x in range(len(lOrg))] , [x[0] for x in lOrg])
    plt.ylim(0, maxprots + maxprots * 0.33)
    plt.title(title)
    plt.legend(loc='best')
    plt.savefig('%s.png'%fname)

    logger.info('%s graph saved (%s.png)'%(title, fname))
    
def plotPanGenome(core, acc, uni):
    plt.clf()
    colors=('r','b','g')
    patches = plt.pie([core, acc, uni], colors=colors,
                      explode=(0.1,0.01,0.01),
                      autopct='%1.1f%%',
                      shadow=True)
    plt.legend((patches[0][0], patches[0][1], patches[0][2]),
               ('Core','Accessory','Unique'),
               loc=(0,-.1))
    plt.title('PanGenome shape')
    plt.savefig('pangenome_shape.png')
    
    logger.info('PanGenome shape graph saved (pangenome_shape.png)')

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