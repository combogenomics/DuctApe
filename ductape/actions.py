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
from ductape.common.utils import slice_it, rgb_to_hex, xstr
from ductape.storage.SQLite.database import DBBase, Project, Genome, Organism, \
    Kegg, Biolog
from matplotlib import cm
import logging
import matplotlib.colors as pltcls
import matplotlib.pyplot as plt
import numpy as np
import os
# No country for warnings
np.seterr(all='ignore')
#

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('ductape.actions')

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

def dAdd(project, orgID, name='', descr='', color=None):
    '''
    Add a single organism
    '''
    org = Organism(project)
    
    # If trying to override a present organism, throw an error
    if org.isOrg(orgID):
        logger.warning('Organism %s is already present'%orgID)
        logger.warning('Remove it before addition')
        return False
    
    if not color:
        org.addOrg(orgID, name=name, description=descr)
    else:
        org.addOrg(orgID, name=name, description=descr, color=color)
    
    logger.info('Added organism %s'%orgID)
    
    return True

def dGenomeAdd(project, orgID, filename):
    '''
    Add a single genome
    '''
    if not os.path.exists(filename):
        logger.error('Fasta file %s may not be present'%(filename))
        return False
    
    filename = os.path.abspath(filename)
    org = Organism(project)
    if not org.isOrg(orgID):
        logger.warning('Organism %s is not present yet!'%orgID)
        return False
    
    gen = Genome(project)
    gen.addProteome(orgID, filename)
    logger.info('Added genome %s, having %d proteins'%
                (orgID, gen.howMany(orgID)))
    return True

def dPhenomeAdd(project, orgID, filename):
    '''
    Add a single phenome
    '''
    from ductape.phenome.biolog import BiologParser, Plate
    
    if not os.path.exists(filename):
        logger.error('Phenomic file %s may not be present'%(filename))
        return False
    
    org = Organism(project)
    if not org.isOrg(orgID):
        logger.warning('Organism %s is not present yet!'%orgID)
        return False
    
    filename = os.path.abspath(filename)
    
    bparser = BiologParser(filename)
    bparser.parse()
    
    if len(bparser.plates) == 0:
        logger.warning('No biolog data was found!')
        return False
    
    # Check the organisms id inside the biolog files
    strainNumbers = set([plate.strainNumber for plate in bparser.plates])
    strainNames = set([plate.strainName for plate in bparser.plates])
    samples = set([plate.sample for plate in bparser.plates])
    
    if orgID not in samples:
        logger.debug('No sign of %s in sample field'%orgID)
    if orgID not in strainNames:
        logger.debug('No sign of %s in strainName field'%orgID)
    if orgID not in strainNumbers:
        logger.debug('No sign of %s in strainNumber field'%orgID)
    
    # TODO: regular expression search
    if orgID in samples:
        if len(samples) > 1:
            logger.warning('''More than one organism ID may be present in this phenomic data file!''')
            logger.warning('''%s'''%' '.join(samples))
            return False
        
        for plate in bparser.plates:
            plate.strain = plate.sample
            
    elif orgID in strainNames:
        if len(strainNames) > 1:
            logger.warning('''More than one organism ID may be present in this phenomic data file!''')
            logger.warning('''%s'''%' '.join(strainNames))
            return False
        
        for plate in bparser.plates:
            plate.strain = plate.strainName
        
    elif orgID in strainNumbers:
        if len(strainNumbers) > 1:
            logger.warning('''More than one organism ID may be present in this phenomic data file!''')
            logger.warning('''%s'''%' '.join(strainNumbers))
            return False
        
        for plate in bparser.plates:
            plate.strain = plate.strainNumber
        
    else:
        logger.warning('''The organism ID you provided was not found inside the phenomic data file''')
        logger.info('''Using it anyway to add this data''')
    
    # Prepare a series of Plate objects to catch the replicas
    # (replicas will be handled by the db tough)
    dPlates={}
    for plate in bparser.plates:
        if plate.plate_id not in dPlates:
            dPlates[plate.plate_id] = Plate(plate.plate_id)
        dPlates[plate.plate_id].addData(plate.strain, plate)
    
    # Grep the wells
    wells = [w for plate in dPlates.itervalues() for w in plate.getWells()]
    
    # Add to the project
    biolog = Biolog(project)
    biolog.addWells(wells, clustered=False)
    
    logger.info('Added phenome %s, having %d biolog plates (%d wells)'%
                (orgID, len(dPlates), len(wells)))
    
    return True

def dPhenomeMultiAdd(project, filename):
    '''
    Add a single phenomic file with multiple organisms in it
    '''
    from ductape.phenome.biolog import BiologParser, Plate
    
    if not os.path.exists(filename):
        logger.error('Phenomic file %s may not be present'%(filename))
        return False
    
    filename = os.path.abspath(filename)
    
    bparser = BiologParser(filename)
    bparser.parse()
    
    if len(bparser.plates) == 0:
        logger.warning('No biolog data was found!')
        return False
    
    # Check the organism ids inside the biolog files
    # Assuming the names are correct AND stored inside the strainName field
    logger.debug('Assuming organism IDs are correct and inside the field strainName')
    strainNames = set([plate.strainName for plate in bparser.plates])
    
    strainNames.discard(None)
    strainNames.discard('')
    
    if len(strainNames) == 0:
        logger.warning('''Field strainName doesn't contain any value (%s)'''%filename)
        return False
        
    logger.info('Found the following organism IDs: %s'%' '.join(strainNames))
    
    for plate in bparser.plates:
        plate.strain = plate.strainName
    
    # TODO: regular expressions verification
    
    orgs = strainNames
    
    for orgID in orgs:
        org = Organism(project)
        if not org.isOrg(orgID):
            logger.warning('Organism %s is not present yet! Skipping...'%orgID)
            continue
        
        # Prepare a series of Plate objects to catch the replicas
        # (replicas will be handled by the db tough)
        dPlates={}
        for plate in bparser.plates:
            if plate.strain == orgID:
                if plate.plate_id not in dPlates:
                    dPlates[plate.plate_id] = Plate(plate.plate_id)
                dPlates[plate.plate_id].addData(plate.strain, plate)
        
        # Grep the wells
        wells = [w for plate in dPlates.itervalues() 
                 for w in plate.getWells()]
        
        # Add to the project
        biolog = Biolog(project)
        biolog.addWells(wells, clustered=False)
        
        logger.info('Added phenome %s, having %d biolog plates (%d wells)'%
                    (orgID, len(dPlates), len(wells)))
    
    return True

def dRemove(project, organisms):
    '''
    Remove all the organism info regarding a particular organism ID(s)
    '''
    org = Organism(project)
    for orgID in organisms:
        if not org.isOrg(orgID):
            logger.warning('Organism %s is not present: skipping'%orgID)
            continue
        
        muts = [mutID for mutID in org.getOrgMutants(orgID)]
        
        org.delOrg(orgID, True)
        logger.info('Successfully removed organism %s'%orgID)
        if len(muts) > 0:
            logger.info('Removed also %d %s mutant(s)'%(len(muts),orgID))
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

def dClear(project, keeporg=False, keepkegg=False):
    '''
    Clear all the organisms data
    '''
    
    if not keeporg:
        org = Organism(project)
        org.delAllOrgs(True)
        logger.info('Successfully removed all organisms data')
    
    if not keepkegg:
        kegg = Kegg(project)
        kegg.clear()
        logger.info('Successfully removed all KEGG data')
    
    return True

def dGenomeClear(project):
    '''
    Clear the genomic tables
    '''
    gen = Genome(project)
    gen.clearAllGenome()
    logger.info('Successfully removed all genomic data')
    return True

def dPhenomeRemove(project, organisms):
    '''
    Remove all the phenomic data about specific organism ID(s)
    '''
    biolog = Biolog(project)
    oCheck = Organism(project)
    for org in organisms:
        if not oCheck.isOrg(org):
            logger.warning('Phenome %s is not present: skipping'%org)
            continue
        biolog.delOrg(org)
        logger.info('Successfully removed phenome %s'%org)
        
    if biolog.atLeastOneParameter():
        logger.warning('The activity must be recalculated')
        
    return True

def dPhenomeClear(project):
    '''
    Clear the phenomic tables
    '''
    biolog = Biolog(project)
    biolog.clearAllPhenome()
    logger.info('Successfully removed all phenomic data')
    return True

def dGenomeDirAdd(project, folder, extension):
    '''
    Add a series of genomes contained in a directory
    '''
    if not os.path.exists(folder):
        logger.error('Fasta folder %s may not be present'%(folder))
        return False
    logger.info('Looking for files with extension %s'%extension)
    
    org = Organism(project)
    
    added = 0
    for infile in os.listdir(folder):
        if infile.split('.')[-1] != extension:
            logger.debug('Skipping file %s'%infile)
            continue
        
        orgID = infile.split('.')[0]
        filename = os.path.join(folder, infile)
        if os.path.isdir(filename):
            continue
        
        if not org.isOrg(orgID):
            logger.warning('Organism %s is not present yet! Skipping...'%orgID)
            continue
        
        if not org.isMutant(orgID):
            if not dGenomeAdd(project, orgID, filename):
                logger.error('Could not add genome %s'%infile)
                return False
        else:
            if not dGenomeMutAdd(project, orgID, filename):
                logger.error('Could not add genome %s'%infile)
                return False
        added += 1
    if added > 0:
        logger.info('Added %d genomes from %s'%
                (added, folder))
    else:
        logger.warning('No genomes were added from %s'%folder)
    return True
    
def dPhenomeDirAdd(project, folder, extension):
    '''
    Add a series of phenomes contained in a directory
    '''
    if not os.path.exists(folder):
        logger.error('Phenomes folder %s may not be present'%(folder))
        return False
    else:
        logger.info('Looking for files with extension %s'%extension)
        
        added = 0
        for infile in os.listdir(folder):
            if infile.split('.')[-1] != extension:
                logger.debug('Skipping file %s'%infile)
                continue
            
            filename = os.path.join(folder, infile)
            if os.path.isdir(filename):
                continue
            
            if dPhenomeMultiAdd(project, filename):
                added += 1
        
        if added > 0:
            logger.info('Added %d phenomic data files from %s'%
                    (added, folder))
        else:
            logger.warning('No phenomes were added from %s'%folder)
        return True

def dMutAdd(project, mutID, mutparent,kind, name='', descr='', color=None):
    '''
    Check and add a mutant
    '''
    org = Organism(project)
    
    if org.isOrg(mutID):
        logger.warning('Organism %s is already present'%mutID)
        logger.warning('Remove it before addition')
        return False
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
    if not color:
        org.addOrg(mutID, name=name, description=descr,
               mutant=True, reference=mutparent, mkind=kind)
    else:
        org.addOrg(mutID, name=name, description=descr,
               mutant=True, reference=mutparent, mkind=kind, color=color)
    
    logger.info('Mutant %s (%s) added'
                %(mutID, org.getOrg(mutID).mkind))
    return True

def dGenomeMutAdd(project, mutID, mutfasta):
    '''
    Check and add a mutant
    '''
    if not os.path.exists(mutfasta):
        logger.error('Fasta file %s may not be present'%(mutfasta))
        return False
    
    org = Organism(project)
    if not org.isOrg(mutID):
        logger.warning('Organism %s is not present yet!'%mutID)
        return False
    
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
            if s[0].lstrip()[0] == '#':
                continue
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

def dPhenomeZero(project, blankfile=None):
    '''
    Takes all the biolog data available and performs the signals zero subtraction
    If blankfile is provided, "blank plates" are parsed and then may be used
    for zero subtraction
    '''
    from ductape.phenome.biolog import BiologParser, getSinglePlates
    from ductape.phenome.biolog import zeroPlates, BiologZero
    
    if blankfile:
        logger.info('Going to use a blank file for zero subtraction')
        
        if not os.path.exists(blankfile):
            logger.error('Blank file %s may not be present'%(blankfile))
            return False
        
        bparser = BiologParser(blankfile)
        bparser.parse()
        
        if len(bparser.plates) == 0:
            logger.warning('The blank file contains no plates!')
            return False
    
    biolog = Biolog(project)
    # Fetch the signals to be subtracted, then convert them
    # to appropriate objects
    sigs = [s for s in biolog.getZeroSubtractableSignals()]
    plates = []
    discarded = set()
    for p in getSinglePlates(sigs):
        if p.plate_id in zeroPlates:
            if blankfile:
                for zp in bparser.plates:
                    if zp.plate_id == p.plate_id:
                        plates.append(p) 
            else:
                plates.append(p)
        else:
            discarded.add(p.plate_id)
    
    if len(plates) == 0:
        logger.warning('No plates can be zero subtracted!')
        logger.warning('Found these plates: %s'%' '.join(discarded))
        return False
    
    if blankfile:
        zsub = BiologZero(plates, blank = True, blankData=bparser.plates)
    else:
        zsub = BiologZero(plates)
        
    if not zsub.zeroSubTract():
        logger.warning('Zero subtraction failed!')
        return False
    
    # Grep the wells
    wells = [w for plate in zsub.plates for w in plate.getWells()]
    
    # Add to the project
    biolog = Biolog(project)
    biolog.addWells(wells, clustered=False, replace=True)
    
    logger.info('Zero subtraction done on %d plates'%len(plates))
    if biolog.atLeastOneParameter():
        logger.warning('The activity must be recalculated')
    
    return True

def dPhenomePurge(project, policy, delta=1, filterplates=[]):
    from ductape.phenome.biolog import getPlates, Experiment
    
    biolog = Biolog(project)
    
    sigs = [s for s in biolog.getAllWells()]
    plates = [p for p in getPlates(sigs, nonmean=True)]
    # The user may want to purge only some plates
    if len(filterplates) > 0:
        for p in plates:
            if p.plate_id not in filterplates:
                while p in plates:
                    plates.remove(p)
    isZero = biolog.atLeastOneZeroSubtracted()

    if len(plates) == 0:
        logger.warning('No phenomic data available, skipping purging')
        return True
    else:
        logger.info('Purging %d phenomic plates'%len(plates))

    exp = Experiment(plates=plates, zero=isZero)
    
    if not exp.purgeReplicas(delta=delta,policy=policy):
        logger.error('Could not purge the phenomic experiment')
        return False

    # Move the discarded wells
    biolog.moveDiscardedWells(exp.discarded)
    logger.info('Purged %d phenomic experiments'%len(exp.discarded))
    
    return True

def dPhenomeRestore(project, plates=[]):
    biolog = Biolog(project)
    
    if not biolog.atLeastOnePurged():
        logger.warning('No phenomic experiment to be restored')
        return True
    
    howmany = biolog.restoreDiscardedWells(plates)
    
    logger.info('Restored %d phenomic experiments'%howmany)
    
    return True

def dGenomeStats(project, svg=False, doPrint=True):
    # Which project are we talking about?
    kind = dSetKind(project)
    
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
            
        plotMapBars(lOrg, 'Single genomes statistics', 'single', svg)
        
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
 
            plotMapBars(lPanGenome, 'PanGenome statistics', 'pangenome_stats',
                        svg)
            plotPanGenome(core, acc, uni, svg)
    
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
                        '%s'%ref_id, svg)
    
    else:
        logger.info('No statistics can be computed at this time')
        return False

    return True

def dPhenomeStats(project, activity=5, delta=3, svg=False, doPrint=True):
    from ductape.phenome.biolog import getPlates, Experiment
    from ductape.phenome.biolog import getOrder
    from itertools import combinations
    
    # Which project are we talking about?
    kind = dSetKind(project)
    
    organism = Organism(project)
    biolog = Biolog(project)
    
    ############################################################################
    # Overall plots
    
    logger.info('Overall plots')
    # Setup an experiment
    sigs = [s for s in biolog.getAllSignals()]
    plates = [p for p in getPlates(sigs)]
    
    isZero = biolog.atLeastOneZeroSubtracted()
    
    category = {}
    categorder = []
    for c in biolog.getPlateCategs():
        categ = c.category.replace(' ','_').replace('&','and')
        
        if categ not in category:
            category[categ] = set()
        category[categ].add(c.plate_id)
        
        if categ not in categorder:
            categorder.append(categ)
    
    exp = Experiment(plates=plates, zero=isZero,
                     category=category, categorder=categorder)
    
    exp.plot(svg=svg)
    
    ############################################################################
    # Activity distribution
    
    logger.info('Activity distributions')
    
    # Fake plot top get the correct bin centers
    y,binEdges=np.histogram(range(10),bins=10)
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    
    if kind == 'single':
        fig = plt.figure()
    else:
        fig = plt.figure(figsize=(12,6))
    
        logger.debug('Overall activity')
        d = biolog.getActivityDistribution()
        x = []
        for k, v in d.iteritems():
            for i in range(v):
                x.append(k)
                
        if len(x) != 0:                   
            ax = fig.add_subplot(1,2,1)
            y,binEdges=np.histogram(x,bins=10)
            ax.plot(bincenters,y,'-o', color='black', linewidth=2)
            
            ax.set_ylim(0,max([x for x in d.itervalues()]) + 20)
            ax.set_xlim(min(bincenters)-1,max(bincenters)+1)
            
            x0,x1 = ax.get_xlim()
            y0,y1 = ax.get_ylim()
            ax.set_aspect((x1-x0)/(y1-y0))
            
            ax.set_xticks(bincenters)
            ax.set_xticklabels([str(x) for x in range(10)])
            
            ax.xaxis.grid(color='gray', linestyle='dashed')
            ax.set_axisbelow(True)
            
            ax.set_xlabel('Activity', size='small')
            ax.set_ylabel('# of wells', size='small')
            ax.set_title('Overall', size='small')
    
    logger.debug('Single organisms activity')
    
    dcolors = getOrganismsColors(project)
    
    if kind == 'single':
        ax = fig.add_subplot(1,1,1)
    else:
        ax = fig.add_subplot(1,2,2)
        
    maxv = []
    for org in organism.getAll():
        d = biolog.getActivityDistributionByOrg(org.org_id)
        x = []
        for k, v in d.iteritems():
            for i in range(v):
                x.append(k)
                
        if len(x) == 0:
            continue
                
        y,binEdges=np.histogram(x,bins=10)
        ax.plot(bincenters,y,'-o', color=dcolors[org.org_id], linewidth=2,
                alpha=0.66, label=org.org_id)
        
        maxv.append(max([x for x in d.itervalues()]))
        
    ax.set_ylim(0,max(maxv) + 20)
    ax.set_xlim(min(bincenters)-1,max(bincenters)+1)
    
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    ax.set_aspect((x1-x0)/(y1-y0))
    
    ax.set_xticks(bincenters)
    ax.set_xticklabels([str(x) for x in range(10)])
    
    ax.xaxis.grid(color='gray', linestyle='dashed')
    ax.set_axisbelow(True)
        
    ax.set_xlabel('Activity', size='small')
    ax.set_ylabel('# of wells', size='small')
    if kind == 'single':
        ax.set_title('Activity distribution', size='large')
    else:
        ax.set_title('Single organisms', size='small')
        fig.suptitle('Activity distribution', size='large')
    
    plt.legend(loc='best',prop={'size':6})
    
    if svg:
        fname = 'Activity.svg'
    else:
        fname = 'Activity.png'
    
    plt.savefig(fname)
    
    logger.info('Saved activity distribution graph (%s)'%fname)
    
    plt.clf()
    
    fig = plt.figure(figsize=(12,6))
    axid = 1
    
    logger.debug('Zero/NonZero distributions')
    
    for bzero in [False,True]:
        d = biolog.getActivityDistributionByZero(bzero)
        x = []
        for k, v in d.iteritems():
            for i in range(v):
                x.append(k)
        
        if len(x) == 0:
            ax = fig.add_subplot(1,2,axid)
            ax.set_xlabel('Activity', size='small')
            ax.set_ylabel('# of wells', size='small')
            if not bzero:
                ax.set_title('Zero subtracted wells', size='small')
            else:
                ax.set_title('NoZero subtracted wells', size='small')
            axid += 1
            continue
            
        ax = fig.add_subplot(1,2,axid)
        
        maxv = []
        for org in organism.getAll():
            d = biolog.getActivityDistributionByZeroAndOrg(org.org_id, bzero)
            x = []
            for k, v in d.iteritems():
                for i in range(v):
                    x.append(k)
            
            if len(x) == 0:
                continue
            
            y,binEdges=np.histogram(x,bins=10)
            ax.plot(bincenters,y,'-o', color=dcolors[org.org_id], linewidth=2,
                    alpha=0.66, label=org.org_id)
            
            maxv.append(max([x for x in d.itervalues()]))
        
        if len(maxv) == 0:
            continue
        
        ax.set_ylim(0,max(maxv) + 20)
        ax.set_xlim(min(bincenters)-1,max(bincenters)+1)
        
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))
        
        ax.xaxis.grid(color='gray', linestyle='dashed')
        ax.set_axisbelow(True)
        
        ax.set_xticks(bincenters)
        ax.set_xticklabels([str(x) for x in range(10)])
            
        ax.set_xlabel('Activity', size='small')
        ax.set_ylabel('# of wells', size='small')
        if not bzero:
            ax.set_title('Zero subtracted wells', size='small')
        else:
            ax.set_title('NoZero subtracted wells', size='small')
        
        axid += 1
        
        plt.legend(loc='best',prop={'size':6})
        
    fig.suptitle('Activity distribution by categories', size='large')
    
    plt.legend(loc='best',prop={'size':6})
    
    if svg:
        fname = 'ActivityZero.svg'
    else:
        fname = 'ActivityZero.png'
        
    plt.savefig(fname)
    
    logger.info('Saved Zero/NoZero activity distribution (%s)'%fname)
    
    plt.clf()
    
    fig = plt.figure(figsize=(24,12))
    axid = 1
    
    logger.debug('Category distributions')
    
    categs = []
    for c in biolog.getPlateCategs():
        if c.category not in categs:
            categs.append(c.category)
    
    for categ in categs:
        d = biolog.getActivityDistributionByCateg(categ)
        x = []
        for k, v in d.iteritems():
            for i in range(v):
                x.append(k)
        
        if len(x) == 0:
            ax = fig.add_subplot(2,4,axid)
            ax.set_xlabel('Activity', size='small')
            ax.set_ylabel('# of wells', size='small')
            ax.set_title('%s'%categ.replace(' ','_').replace('&','and'), size='small')
            axid += 1
            continue
            
        ax = fig.add_subplot(2,4,axid)
        
        logger.debug('Plotting category %s activities'%categ)
        maxv = []
        for org in organism.getAll():
            d = biolog.getActivityDistributionByCategAndOrg(categ, org.org_id)
            x = []
            for k, v in d.iteritems():
                for i in range(v):
                    x.append(k)
            
            if len(x) == 0:
                continue
            
            y,binEdges=np.histogram(x,bins=10)
            ax.plot(bincenters,y,'-o', color=dcolors[org.org_id], linewidth=2,
                    alpha=0.66, label=org.org_id)
            
            maxv.append(max([x for x in d.itervalues()]))
        
        if len(maxv) == 0:
            continue
        
        ax.set_ylim(0,max(maxv) + 20)
        ax.set_xlim(min(bincenters)-1,max(bincenters)+1)
        
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))
        
        ax.xaxis.grid(color='gray', linestyle='dashed')
        ax.set_axisbelow(True)
        
        ax.set_xticks(bincenters)
        ax.set_xticklabels([str(x) for x in range(10)])
            
        ax.set_xlabel('Activity', size='small')
        ax.set_ylabel('# of wells', size='small')
        ax.set_title('%s'%categ.replace(' ','_').replace('&','and'), size='small')
        
        axid += 1
        
        plt.legend(loc='best',prop={'size':6})
        
    fig.suptitle('Activity distribution by categories', size='large')
    
    plt.legend(loc='best',prop={'size':6})
    
    if svg:
        fname = 'ActivityCateg.svg'
    else:
        fname = 'ActivityCateg.png'
    
    plt.savefig(fname)
    
    logger.info('Saved category activity distribution (%s)'%fname)
    
    plt.clf()
    
    ############################################################################
    # Activity Boxplots
    logger.info('Activity boxplots')
    
    if kind == 'single':
        fig = plt.figure()
    else:
        fig = plt.figure(figsize=(12,6))
    
        d = biolog.getActivityDistribution()
        x = []
        for k, v in d.iteritems():
            for i in range(v):
                x.append(k)
                
        if len(x) != 0:                   
            ax = fig.add_subplot(1,2,1)
            bplot=ax.boxplot(x,1,vert=0)
            
            colorBoxPlot(ax, bplot, ['black'])
            
            ax.set_xlim(-1,10)
            ax.set_xticks(range(10))
            
            ax.get_yaxis().set_ticks([])
            
            x0,x1 = ax.get_xlim()
            y0,y1 = ax.get_ylim()
            ax.set_aspect((x1-x0)/(y1-y0))
            
            ax.xaxis.grid(color='gray', linestyle='dashed')
            ax.set_axisbelow(True)
            
            ax.set_xlabel('Activity', size='small')
            ax.set_title('Overall', size='small')
    
    logger.debug('Single organisms activity')
    
    dcolors = getOrganismsColors(project)
    
    if kind == 'single':
        ax = fig.add_subplot(1,1,1)
    else:
        ax = fig.add_subplot(1,2,2)
        
    # Get the organisms order (to have a nice order in case of mutants)
    orgs = []
    # Orgs actually used
    borgs = []
    bcolors = [] 
    data = []
    if kind == 'mutants':
        refs = [org.org_id
                    for org in organism.getAll()
                    if not organism.isMutant(org.org_id)]
        
        for ref_id in refs:
            orgs.append(ref_id)            
            for x in organism.getOrgMutants(ref_id):
                orgs.append(x)
    else:
        for org in organism.getAll():
            orgs.append(org.org_id)
    
    for org_id in orgs:
        d = biolog.getActivityDistributionByOrg(org_id)
        x = []
        for k, v in d.iteritems():
            for i in range(v):
                x.append(k)
                
        if len(x) == 0:
            continue
        borgs.append(org_id)
        bcolors.append(dcolors[org_id])
        data.append(x)
    
    borgs = borgs[::-1]
    bcolors = bcolors[::-1]
    data = data[::-1]
    
    bplot=ax.boxplot(data,1,vert=0)
            
    colorBoxPlot(ax, bplot, bcolors)
    
    if len(borgs) > 0:
        ax.set_xlim(-1,10)
        ax.set_xticks(range(10))
        
        ax.set_yticklabels(borgs, size='x-small')
        
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))
        
        ax.xaxis.grid(color='gray', linestyle='dashed')
        ax.set_axisbelow(True)
        
        ax.set_xlabel('Activity', size='small')
        ax.set_ylabel('Organisms', size='x-small')
    
        if kind == 'single':
            ax.set_title('Activity boxplots', size='large')
        else:
            ax.set_title('Single organisms', size='small')
            fig.suptitle('Activity boxplots', size='large')
    
        if svg:
            fname = 'ActivityBoxplot.svg'
        else:
            fname = 'ActivityBoxplot.png'
            
        plt.savefig(fname)
            
        logger.info('Saved activity boxplots (%s)'%fname)
    
    plt.clf()
    
    fig = plt.figure(figsize=(24,12))
    axid = 1
    
    logger.debug('Category distributions')
    
    categs = []
    for c in biolog.getPlateCategs():
        if c.category not in categs:
            categs.append(c.category)
    
    for categ in categs:
        d = biolog.getActivityDistributionByCateg(categ)
        x = []
        for k, v in d.iteritems():
            for i in range(v):
                x.append(k)
        
        if len(x) == 0:
            ax = fig.add_subplot(2,4,axid)
            ax.set_xlabel('Activity', size='small')
            ax.set_ylabel('Organisms', size='x-small')
            ax.set_title('%s'%categ.replace(' ','_').replace('&','and'), size='small')
            axid += 1
            continue
            
        ax = fig.add_subplot(2,4,axid)
        
        logger.debug('Plotting category %s activities'%categ)
        
        # Orgs actually used
        borgs = []
        bcolors = [] 
        data = []
        for org_id in orgs:
            d = biolog.getActivityDistributionByCategAndOrg(categ, org_id)
            x = []
            for k, v in d.iteritems():
                for i in range(v):
                    x.append(k)
                    
            if len(x) == 0:
                continue
            borgs.append(org_id)
            bcolors.append(dcolors[org_id])
            data.append(x)
        
        borgs = borgs[::-1]
        bcolors = bcolors[::-1]
        data = data[::-1]
        
        bplot=ax.boxplot(data,1,vert=0)
                
        colorBoxPlot(ax, bplot, bcolors)
        
        if len(borgs) > 0:
            ax.set_xlim(-1,10)
            ax.set_xticks(range(10))
            
            ax.set_yticklabels(borgs, size='x-small')
            
            x0,x1 = ax.get_xlim()
            y0,y1 = ax.get_ylim()
            ax.set_aspect((x1-x0)/(y1-y0))
            
            ax.xaxis.grid(color='gray', linestyle='dashed')
            ax.set_axisbelow(True)
            
            ax.set_xlabel('Activity', size='small')
            ax.set_ylabel('Organisms', size='x-small')
            ax.set_title('%s'%categ, size='small')
        
        axid += 1
            
    fig.suptitle('Activity boxplots by categories', size='large')
    
    if svg:
        fname = 'ActivityCategBoxplot.svg'
    else:
        fname = 'ActivityCategBoxplot.png'
    
    plt.savefig(fname)
    
    logger.info('Saved category activity boxplots (%s)'%fname)
    
    plt.clf()
    
    ############################################################################
    # Statistics printing
    logger.info('Active wells stats')
    
    # Get the organisms order (to have a nice order in case of mutants)
    orgs = []
    refs = {}
    if kind == 'mutants':
        refIDs = [org.org_id
                    for org in organism.getAll()
                    if not organism.isMutant(org.org_id)]
        
        for ref_id in refIDs:
            refs[ref_id] = []
            orgs.append(ref_id)            
            for x in organism.getOrgMutants(ref_id):
                orgs.append(x)
                refs[ref_id].append(x)
    else:
        for org in organism.getAll():
            orgs.append(org.org_id)
            
    header = 'Active wells (%% of wells with activity >= %d)'%activity
    if doPrint:
        print header
    else:
        logger.info(header)
    
    header = '\t'.join( ['Category'] + orgs )
    if doPrint:
        print header
    else:
        logger.info(header)
        
    for categ in categorder:
        line = [categ]
        for org_id in orgs:
            wells = filter(lambda x:x.plate_id in category[categ],
                                [w for w in exp.getAverageWells(org_id)])
            total = float(len(wells))
            active = float(len(filter(lambda x: x.activity >= activity,
                                      wells)))
            
            try:
                line.append( str( (active/total) * 100) )
            except:
                line.append('N/A')
    
        line = '\t'.join(line)
        if doPrint:
            print line
        else:
            logger.info(line)
    
    if kind == 'single':
        return True
    
    logger.info('Active differences stats')
    
    header = '\t'.join( ['Category', 'Average difference',
             'Main differences (%% of wells whose average difference >= %d)'%delta] )
    if doPrint:
        print header
    else:
        logger.info(header)
        
    for categ in categorder:
        pwlist = filter(lambda x: x[0] in category[categ],
                    [pw for pw in getOrder()])
        
        pwdiff = []
        pwavgdiff = []
        pwcount = 0
        if kind == 'mutants':
            for pid, wid in pwlist:
                try:
                    # For simplicity reasons, we use the overall differences
                    diff = []
                    for ref_id in refIDs:
                        for mut_id in refs[ref_id]:
                            if ref_id not in exp.sumexp[pid][wid]:continue
                            if mut_id not in exp.sumexp[pid][wid]:continue
                            
                            pwdiff.append( exp.sumexp[pid][wid][mut_id].activity -
                                           exp.sumexp[pid][wid][ref_id].activity )
                            diff.append( abs(
                                       exp.sumexp[pid][wid][mut_id].activity -
                                       exp.sumexp[pid][wid][ref_id].activity) )
                            
                    pwavgdiff.append(np.array(diff).mean())
                    pwcount += 1
                except:pass
        else:
            for pid, wid in pwlist:
                try:
                    diff = []
                    for oid, oid1 in combinations(orgs, 2):
                        if oid not in exp.sumexp[pid][wid]:continue
                        if oid1 not in exp.sumexp[pid][wid]:continue
                                               
                        pwdiff.append( abs(
                                   exp.sumexp[pid][wid][oid].activity -
                                   exp.sumexp[pid][wid][oid1].activity) )
                        diff.append( abs(
                                   exp.sumexp[pid][wid][oid].activity -
                                   exp.sumexp[pid][wid][oid1].activity) )
                    
                    pwavgdiff.append(np.array(diff).mean())
                    pwcount += 1 
                except:pass
        
        total = float(pwcount)
        if len(pwdiff) == 0:
            avgdiff = 'N/A'
        else:
            avgdiff = str(np.array(pwdiff).mean())
        overdiff = float(len(filter(lambda x: x >= delta, pwavgdiff)))
        try:
            maindiff = str((overdiff / total) * 100)
        except:
            maindiff = 'N/A'
        
        line = '\t'.join( [categ, avgdiff, maindiff] )
        if doPrint:
            print line
        else:
            logger.info(line)
       
    ############################################################################
    
    if kind == 'single':
        return True
    
    logger.info('Unique metabolic functions')
    
    header = 'Unique metabolic functions (%% of wells with delta >= %d)'%delta
    if doPrint:
        print header
    else:
        logger.info(header)
    
    header = '\t'.join( [''] + orgs )
    if doPrint:
        print header
    else:
        logger.info(header)
    
    UpUnique = {}
    for oid in orgs:
        UpUnique[oid] = 0
    DownUnique = {}
    for oid in orgs:
        DownUnique[oid] = 0
    pwcount = 0
    for pid, wid in getOrder():
        try:
            # Sort by activity
            acts = sorted([exp.sumexp[pid][wid][oid] for oid in orgs],
                          key=lambda x: x.activity)
            
            amax = acts[-1]
            amaxindex = acts.index(amax)
            if amax.activity - acts[amaxindex - 1].activity >= delta:
                UpUnique[amax.strain] += 1

            amin = acts[0]
            if acts[1].activity - amin.activity >= delta:
                DownUnique[amin.strain] += 1
            
            pwcount += 1
        except:pass
        
    total = float(pwcount)
    
    line = '\t'.join(['More active'] + [str(UpUnique[oid]) for oid in orgs])    
    if doPrint:
        print line
    else:
        logger.info(line)
    
    line = '\t'.join(['Less active'] + [str(DownUnique[oid]) for oid in orgs])    
    if doPrint:
        print line
    else:
        logger.info(line)
            
    line = '\t'.join(['Total'] + [str(UpUnique[oid] + DownUnique[oid]) for oid in orgs])    
    if doPrint:
        print line
    else:
        logger.info(line)
        
    line = '\t'.join(['%'] + [str(((UpUnique[oid] + DownUnique[oid])/total)*100) for oid in orgs])    
    if doPrint:
        print line
    else:
        logger.info(line)
    
    return True

def dPhenomeRings(project, delta=1, difforg=None, svg=False):
    from ductape.phenome.biolog import getPlates, Experiment
    from ductape.phenome.biolog import getOrder
    
    # Which project are we talking about?
    kind = dSetKind(project)
    
    if kind == 'mutants' and difforg:
        logger.warning('Reference organism(s) will be used for diff mode')
        difforg = None
    
    organism = Organism(project)
    if difforg:
        if not organism.isOrg(difforg):
            logger.error('Organism %s is not present yet!'%difforg)
            return False
        logger.info('Diff mode: using organism %s as reference'%difforg)
    
    biolog = Biolog(project)
    
    # Setup an experiment
    sigs = [s for s in biolog.getAllSignals()]
    plates = [p for p in getPlates(sigs)]
    
    isZero = biolog.atLeastOneZeroSubtracted()
    
    category = {}
    categorder = []
    for c in biolog.getPlateCategs():
        categ = c.category.replace(' ','_').replace('&','and')
        
        if categ not in category:
            category[categ] = set()
        category[categ].add(c.plate_id)
        
        if categ not in categorder:
            categorder.append(categ)
    
    exp = Experiment(plates=plates, zero=isZero,
                     category=category, categorder=categorder)
    
    ############################################################################
    # Activity rings (!!!)
    # Thanks to stackoverflow for that
    # http://stackoverflow.com/questions/12803883
    
    logger.info('Activity ring')
    
    fig = plt.figure(figsize=(25,25), dpi=300)
    # Polar plot!
    ax = fig.add_subplot(111, polar = True)
    # Start from the top
    ax.set_theta_offset(np.pi/2)
    # Go clockwise
    ax.set_theta_direction(-1)
    
    # Get the organisms order (to have a nice order in case of mutants)
    orgs = []
    muts = {}
    if kind == 'mutants':
        refs = [org.org_id
                    for org in organism.getAll()
                    if not organism.isMutant(org.org_id)]
        
        for ref_id in refs:
            orgs.append(ref_id)            
            for x in organism.getOrgMutants(ref_id):
                orgs.append(x)
                muts[x] = ref_id
    elif difforg:
        for org in organism.getAll():
            orgs.append(org.org_id)
            if org.org_id != difforg:
                muts[org.org_id] = difforg
    else:
        for org in organism.getAll():
            orgs.append(org.org_id)
    
    # "Legend"
    i = 0.1
    for org_id in orgs:
        radius = np.linspace(i, i+0.1, 10)
        theta = np.linspace(0, 2*np.pi, 628)
        R,T  = np.meshgrid(radius,theta)
    
        ax.pcolor(T, R, [[1 for y in range(10)] for x in range(628)], cmap=cm.Greys,
                  vmin=0, vmax=9)
        
        ax.text(0, i+0.03, org_id, size=17, weight='black', alpha=0.77, ha='center')
        
        i += 0.15
        
    # Category and plate/well fail-safe tweaks
    # What if a plate/well is missing?
    # What is the category order?
    categpworder = {}
    for categ in categorder:
        categpworder[categ] = []
        for pid, wid in getOrder(sorted(category[categ])):
            # Check what we can discard (plates)
            if pid in exp.sumexp:
                categpworder[categ].append((pid, wid))
    
    # Categ colors
    from itertools import cycle
    categcolor = {}
    colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
    for categ, col in zip(categorder, colors):
        categcolor[categ] = col
    
    # Starting point 
    if i < 0.5:
        i = 0.5
    for org_id in orgs:
        acts = []
        
        for categ in categorder:
            if len(categpworder[categ]) == 0:continue
            
            for pid, wid in categpworder[categ]:
                try:
                    w = exp.sumexp[pid][wid][org_id]
                except:
                    acts.append(np.nan)
                    continue
                if (kind == 'mutants' or difforg) and org_id in muts:
                    # Check if the reference has a value for this
                    try:
                        ref_id = muts[org_id]
                        refact = exp.sumexp[w.plate_id][w.well_id][ref_id].activity
                        if abs(w.activity - refact) <= delta:
                            acts.append(np.nan)
                        else:
                            acts.append(w.activity - refact)
                    except Exception, e:
                        logger.warning(e)
                        acts.append(np.nan)
                else:
                    if w.activity is None:
                        acts.append(np.nan)
                    else:
                        acts.append(w.activity)
                    
        radius = np.linspace(i, i+0.2, 10)
        theta = np.linspace(0, 2*np.pi, len(acts))
        R,T  = np.meshgrid(radius,theta)
        
        if (kind == 'mutants' or difforg) and org_id in muts:
            cmap = cm.PuOr
            vmin = -9
            vmax = 9
        else:
            cmap = cm.RdYlGn
            vmin = 0
            vmax = 9
            
        cmap.set_under('#F8F8F8',1.)
        
        ax.pcolor(T, R, [[x for y in range(10)] for x in acts], cmap=cmap,
                  vmin=vmin, vmax=vmax)
        
        i += 0.25

    # Categ archs
    # Total points
    total = 0
    for pw in categpworder.itervalues():
        total += len(pw)
    total = float(total)
    
    start_arch = 0
    stop_arch = 0
    for categ in categorder:
        if len(categpworder[categ]) == 0:continue
        
        # This category proportion over the others
        categprop = float(len(categpworder[categ])) / total
        
        stop_arch += 2*np.pi*categprop
        theta = np.linspace(start_arch, stop_arch, 100)
        
        ax.plot(theta, [i for t in theta], color=categcolor[categ],
                linewidth=20, label=categ)
        
        start_arch += 2*np.pi*categprop

    #Turn off polar labels
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    ax.set_rmax(i)

    # Colorbar
    import matplotlib.colors as colors
    cNorm  = colors.Normalize(vmin=0, vmax=9)
    scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cm.RdYlGn)
    scalarMap.set_array(np.array(range(10)))
    cax = fig.add_axes([0.93, 0.2, 0.03, 0.6])
    cax.text(0.50, 1.01, 'Activity', size=20, ha='center')
    plt.colorbar(scalarMap, cax=cax)
    
    if (kind == 'mutants' or difforg):
        cNorm  = colors.Normalize(vmin=-9, vmax=9)
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cm.PuOr)
        scalarMap.set_array(np.array(range(-9,9,20)))
        cax = fig.add_axes([0.04, 0.2, 0.03, 0.6])
        cax.text(0.50, 1.01, 'Delta activity', size=20, ha='center')
        plt.colorbar(scalarMap, cax=cax)
        
        # Title
        ax.set_title('Activity ring (diff mode)', size=35)
    else:
        # Title
        ax.set_title('Activity ring', size=35)
    
    ax.legend(loc='best')
    
    if svg:
        fname = 'ActivityRing.svg'
    else:
        fname = 'ActivityRing.png'
    
    plt.savefig(fname)
    
    logger.info('Saved activity ring (%s)'%fname)
    
    plt.clf()
    
    return True

def dKeggImport(project, infile):
    kegg = Kegg(project)
    
    logger.info('Importing KEGG metabolic network')
    
    kegg.importKegg(open(infile))
    
    i = 0
    for l in open(infile):
        i += 1
    logger.info('Imported %d KEGG metabolic network entries'%i)
    
    return True

def dKeggExport(project):
    kegg = Kegg(project)
    
    logger.info('Exporting KEGG metabolic network')
    
    fname = 'kegg.tsv'
    fout = open(fname,'w')
    
    i = 0
    for row in kegg.exportKegg():
        fout.write(row + '\n')
        i += 1
        
    fout.close()
    
    logger.info('Exported %d KEGG entries (%s)'%(i, fname))
    
    return True

def dGenomeExport(project):
    from Bio import SeqIO
    
    # Is there something to be exported?
    organism = Organism(project)
    
    if organism.howMany() == 0:
        logger.info('No genomic data can be exported at this time')
        return False
    
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
        fout.write('#%s\t%s\n'%('prot_id', 'ko_id'))
        i = 0
        for prot_id, ko_id in kegg.getAllKO(org.org_id):
            fout.write('%s\t%s\n'%(prot_id, ko_id.lstrip('ko:')))
            i += 1
        fout.close()
        
        if i == 0:
            os.remove(fname)
            logger.warning('No KO links available for %s'%org.org_id)
        else:
            logger.info('Saved %d KO links for %s (%s)'%(i, org.org_id,
                                                     fname))
        
    logger.info('Exporting Kegg reactions data')
    
    for org in organism.getAll():
        fname = 'reactions_%s.tsv'%org.org_id
        fout = open(fname,'w')
        fout.write('#%s\t%s\n'%('prot_id', 're_id'))
        i = 0
        for prot_id, re_id in kegg.getAllReactions(org.org_id):
            fout.write('%s\t%s\n'%(prot_id, re_id.lstrip('rn:')))
            i += 1
        fout.close()
        
        if i == 0:
            os.remove(fname)
            logger.warning('No Kegg reactions available for %s'%org.org_id)
        else:
            logger.info('Saved %d Kegg reactions links for %s (%s)'%
                    (i, org.org_id, fname))
        
    proj = Project(project)
    
    if proj.isPanGenome():
        logger.info('Exporting pangenome data')
        
        dG = genome.getPanGenome()
        if len(dG) == 0:
            logger.warning('No pangenome available')
        else:
            fname = 'pangenome.tsv'
            fout = open(fname,'w')
            fout.write('#%s\t%s\n'%('orth_id', 'prot_id'))
            for group, prots in dG.iteritems():
                for prot in prots:
                    fout.write('%s\t%s\n'%(group,prot))
            fout.close()
            
            logger.info('Exported %d orthologs (%s)'%(len(dG),fname))
            
            fname = 'pangenome_category.tsv'
            fout = open(fname,'w')
            fout.write('#%s\t%s\t%s\n'%('orth_id', 'category', 'organism(s)'))
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

def dPhenomeExport(project):
    # Is there something to be exported?
    organism = Organism(project)
    
    if organism.howMany() == 0:
        logger.info('No phenomic data can be exported at this time')
        return False
    
    biolog = Biolog(project)    
    
    # Which project are we talking about?
    kind = dSetKind(project)    
    
    logger.info('Exporting single organism(s) phenomic data')
    
    for org in organism.getAll():
        fname = 'phenome_%s.tsv'%org.org_id
        fout = open(fname,'w')
        fout.write('#' + '\t'.join(['plate_id', 'well_id', 'chemical',
                                'category',
                                'moa', 'co_id', 'replica', 'activity',
                                'min', 'max', 'height', 'plateau', 'slope',
                                'lag', 'area']) + '\n')
        i = 0
        for w in biolog.getOrgWells(org.org_id):
            wdet = biolog.getWell(w.plate_id, w.well_id)
            fout.write('\t'.join([xstr(x) for x in [w.plate_id, w.well_id, wdet.chemical,
                                  wdet.category, wdet.moa, wdet.co_id]] +
                                  [xstr(x) for x in [w.replica, w.activity,
                                                    w.min, w.max, w.height,
                                                    w.plateau, w.slope,
                                                    w.lag, w.area]])
                       + '\n')
            i += 1
        fout.close()
        
        if i == 0:
            os.remove(fname)
            logger.warning('No phenomic experiments available for %s'%org.org_id)
        else:            
            logger.info('Saved %d phenomic experiments from %s (%s)'%(i,
                                                org.org_id,
                                                'phenome_%s.tsv'%org.org_id))
        
        # Export the average activity if we have replicas
        if biolog.howManyReplicasByOrg(org.org_id) > 1:
            fname = 'phenome_avg_%s.tsv'%org.org_id
            fout = open(fname,'w')
            fout.write('#' + '\t'.join(['plate_id', 'well_id', 'chemical',
                                    'category',
                                    'moa', 'co_id', 'replica',
                                    'avg activity']) + '\n')
            i = 0
            for w in biolog.getOrgDistinctWells(org.org_id):
                wdet = biolog.getWell(w.plate_id, w.well_id)
                fout.write('\t'.join([xstr(x) for x in [w.plate_id, w.well_id, wdet.chemical,
                                      wdet.category, wdet.moa, wdet.co_id,
                                      xstr(biolog.getAvgActivity(w.plate_id,
                                                            w.well_id,
                                                            org.org_id))]])
                           + '\n')
                i += 1
            fout.close()            

            if i == 0:
                os.remove(fname)
                logger.warning('No average phenomic experiments available'%org.org_id)
            else:            
                logger.info('Saved %d average phenomic experiments from %s (%s)'%(i,
                                                    org.org_id,
                                                    'phenome_avg_%s.tsv'%org.org_id))
    
    if kind == 'pangenome':
        logger.info('Exporting combined phenomes')
        
        fname = 'phenome_combined.tsv'
        fout = open(fname,'w')
        fout.write('#' + '\t'.join(['', '', '','','', '', '',
                                            'activity']) + '\n')            
        fout.write('#' + '\t'.join(['plate_id', 'well_id', 'chemical',
                                'category',
                                'moa', 'co_id', 'replica'] +
                                [x.org_id for x in organism.getAll()])
                   + '\n')        
        
        i = 0
        for w in biolog.getDistinctWells(replica=True):
            wdet = biolog.getWell(w.plate_id, w.well_id)
            fout.write('\t'.join([xstr(x) for x in [w.plate_id, w.well_id, wdet.chemical,
                                wdet.category, wdet.moa, wdet.co_id,
                                w.replica]] +
                                [xstr(biolog.getOneWell(w.plate_id, w.well_id,
                                                   x.org_id, w.replica).activity)
                                 for x in organism.getAll()]) + '\n')
            i += 1
        fout.close()            

        if i == 0:
            os.remove(fname)
            logger.warning('No combined phenomic experiments available')
        else:            
            logger.info('Saved %d combined phenomic experiments (%s)'%(i,
                        'phenome_combined.tsv'))
            
        # Export the average activity if we have replicas
        if biolog.howManyReplicas() > 1:
            i = 0
            fname = 'phenome_avg_combined.tsv'
            fout = open(fname,'w')
            fout.write('#' + '\t'.join(['', '', '','','', '', '',
                                                'avg activity']) + '\n')            
            fout.write('#' + '\t'.join(['plate_id', 'well_id', 'chemical',
                                    'category',
                                    'moa', 'co_id'] +
                                    [x.org_id for x in organism.getAll()])
                       + '\n')        
            
            for w in biolog.getDistinctWells(replica=False):
                wdet = biolog.getWell(w.plate_id, w.well_id)
                fout.write('\t'.join([xstr(x) for x in [w.plate_id, w.well_id, wdet.chemical,
                                    wdet.category, wdet.moa, wdet.co_id]] +
                                    [xstr(biolog.getAvgActivity(w.plate_id,
                                                                w.well_id,
                                                                x.org_id))
                                     for x in organism.getAll()]) + '\n')
                i += 1
            fout.close()            
    
            if i == 0:
                os.remove(fname)
                logger.warning('No combined average phenomic experiments available')
            else:            
                logger.info('Saved %d combined average phenomic experiments (%s)'%(i,
                            'phenome_avg_combined.tsv'))
    
    elif kind == 'mutants':
        logger.info('Exporting combined phenomes')
        
        refs = [org.org_id
                for org in organism.getAll()
                if not organism.isMutant(org.org_id)]
        
        fname = 'phenome_combined.tsv'
        fout = open(fname,'w')
        fout.write('#' + '\t'.join(['', '', '','','', '', '',
                                            'activity']) + '\n')            
        fout.write('#' + '\t'.join(['plate_id', 'well_id', 'chemical',
                                'category',
                                'moa', 'co_id', 'replica']))

        for ref in refs:
            fout.write('\t' + '\t'.join([ref] + [x for x in organism.getOrgMutants(ref)]))
        fout.write('\n')
        
        i = 0
        for w in biolog.getDistinctWells(replica=True):
            wdet = biolog.getWell(w.plate_id, w.well_id)
            fout.write('\t'.join([xstr(x) for x in [w.plate_id, w.well_id, wdet.chemical,
                                wdet.category, wdet.moa, wdet.co_id,
                                w.replica]]))
            for ref in refs:
                fout.write('\t' + '\t'.join([xstr(biolog.getOneWell(w.plate_id,
                                                w.well_id,
                                                ref,
                                                w.replica).activity)] + 
                                            [xstr(biolog.getOneWell(w.plate_id,
                                                w.well_id,
                                                x,
                                                w.replica).activity)
                                             for x in organism.getOrgMutants(ref)]))
            fout.write('\n')
            i += 1
        fout.close()            

        if i == 0:
            os.remove(fname)
            logger.warning('No combined phenomic experiments available')
        else:            
            logger.info('Saved %d combined phenomic experiments (%s)'%(i,
                        'phenome_combined.tsv'))
        
        # Export the average activity if we have replicas
        if biolog.howManyReplicas() > 1:
            i = 0
            fname = 'phenome_avg_combined.tsv'
            fout = open(fname,'w')
            fout.write('#' + '\t'.join(['', '', '','','', '', '',
                                                'avg activity and deltas']) + '\n')            
            fout.write('#' + '\t'.join(['plate_id', 'well_id', 'chemical',
                                    'category',
                                    'moa', 'co_id']))
    
            for ref in refs:
                fout.write('\t' + '\t'.join([ref] + [x for x in organism.getOrgMutants(ref)]))
            fout.write('\n')
            
            for w in biolog.getDistinctWells(replica=False):
                wdet = biolog.getWell(w.plate_id, w.well_id)
                fout.write('\t'.join([xstr(x) for x in [w.plate_id, w.well_id, wdet.chemical,
                                    wdet.category, wdet.moa, wdet.co_id]]))
                for ref in refs:
                    ref_act = biolog.getAvgActivity(w.plate_id, w.well_id, ref)
                    fout.write('\t' + '\t'.join([xstr(ref_act)] + 
                                                [xstr(ref_act - biolog.getAvgActivity(
                                                    w.plate_id,
                                                    w.well_id,
                                                    x))
                                                 for x in organism.getOrgMutants(ref)]))
                fout.write('\n')
                i += 1
            fout.close()            
    
            if i == 0:
                os.remove(fname)
                logger.warning('No combined average phenomic experiments available')
            else:            
                logger.info('Saved %d combined average phenomic experiments (%s)'%(i,
                            'phenome_avg_combined.tsv'))        
                       
    return True
            
def dSetKind(project):
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
    kind = dSetKind(project)
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

def getPathsComps(project):
    kegg = Kegg(project)
    # Get the pathway - compounds links
    paths = {}
    for pR in kegg.getPathComps():
        if pR.path_id in ['path:rn01100','path:rn01110','path:rn01120']:
            continue
        if pR.path_id not in paths:
            paths[pR.path_id] = []
        paths[pR.path_id].append(pR.co_id)
        
    return paths

def getOrganismsColors(project):
    '''
    Check the colors assigned to the organisms and return a dictionary
    If no colors are assigned, they are assigned automatically
    '''
    organism = Organism(project)
    
    colors = {}
    for org in organism.getAll():
        if not org.color or org.color == '':
            colors[org.org_id] = None
        else:
            colors[org.org_id] = org.color

    orgs = colors.keys()
    for org, color in colors.iteritems():
        # Automatic assignment, probably not the best choiche
        # if we got some organism assigned and some others not
        if not color:
            if len(orgs) == 1:
                autocolor = plt.get_cmap('jet')(float( orgs.index(org) )/(len(orgs)))
            else:
                autocolor = plt.get_cmap('jet')(float( orgs.index(org) )/(len(orgs)-1))
            autocolor = pltcls.rgb2hex(autocolor)
            colors[org] = autocolor
            organism.setColor(org, autocolor)
            logger.info('Automatically assigned color to %s'%org)
    
    return colors

def prepareColors(maximum, colorrange):
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

def createLegend(kind):
    '''
    Create a color scheme legend
    '''
    # TODO: a more centralized color scheme
    fig = plt.figure()
    fname = 'legend.png'
    rmatrix = np.outer(np.arange(0.33,1,0.01),np.ones(7))
    if kind == 'mutants' or 'singlediff':
        cmatrix = np.outer(np.arange(-9,9,0.1),np.ones(7))
    elif kind == 'pangenome':
        cmatrix = np.outer(np.arange(0.33,1,0.01),np.ones(7))
    else:
        cmatrix = np.outer(np.arange(0,9,0.1),np.ones(7))
         
    if kind == 'pangenome':
        ax = fig.add_subplot(141)
        ax.imshow(rmatrix, cmap=cm.Blues, vmin=0, vmax=1)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Core')
        
        ax = fig.add_subplot(142)
        ax.imshow(rmatrix, cmap=cm.Greens, vmin=0, vmax=1)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Core and Disp.')
        
        ax = fig.add_subplot(143)
        ax.imshow(rmatrix, cmap=cm.Oranges, vmin=0, vmax=1)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Dispensable')
        
        ax = fig.add_subplot(144)
        ax.imshow(cmatrix, cmap=cm.Purples, vmin=0, vmax=1)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Mean activity diff.')
        
        fig.savefig(fname)
            
    elif kind == 'single':
        ax = fig.add_subplot(121)
        ax.imshow(rmatrix, cmap=cm.Greens, vmin=0, vmax=1)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Reactions')
        
        ax = fig.add_subplot(122)
        ax.imshow(cmatrix, cmap=cm.RdYlGn, vmin=0, vmax=9)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Phenomic activity')
        
        fig.savefig(fname)
        
    elif kind == 'singlediff':
        ax = fig.add_subplot(121)
        ax.imshow(rmatrix, cmap=cm.Greens, vmin=0, vmax=1)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Reactions')
        
        ax = fig.add_subplot(122)
        ax.imshow(cmatrix, cmap=cm.PuOr, vmin=-9, vmax=9)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Delta-activity')
        
        fig.savefig(fname)

    elif kind == 'mutants':
        ax = fig.add_subplot(141)
        ax.imshow(rmatrix, cmap=cm.Greens, vmin=0, vmax=1)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Wild-type')
        
        ax = fig.add_subplot(142)
        ax.imshow(rmatrix, cmap=cm.copper_r, vmin=0, vmax=1)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Wild-type and Mutated')
        
        ax = fig.add_subplot(143)
        ax.imshow(rmatrix, cmap=cm.Reds, vmin=0, vmax=1)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Mutated')
        
        ax = fig.add_subplot(144)
        ax.imshow(cmatrix, cmap=cm.PuOr, vmin=-9, vmax=9)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Delta-activity')
        
        fig.savefig(fname)
        
    return fname

def colorBoxPlot(ax, bplot, colors):
    '''
    color the boxplots with the desired color
    
    Borrowed from the amazing matplotlib gallery
    http://matplotlib.org/examples/pylab_examples/boxplot_demo2.html
    '''
    import matplotlib.artist as art
    from matplotlib.patches import Polygon
    
    # Caps, Fliers and Whiskers
    tocolor = []
    for i in range(1,len(bplot['whiskers']),2):
        color = colors[i/2]
    
        tocolor.append((bplot['caps'][i-1], color))
        tocolor.append((bplot['caps'][i], color))
        
        tocolor.append((bplot['fliers'][i-1], color))
        tocolor.append((bplot['fliers'][i], color))
        
        tocolor.append((bplot['whiskers'][i-1], color))
        tocolor.append((bplot['whiskers'][i], color))
    
    for obj, color in tocolor:
        art.setp(obj, color=color, linewidth=2)
    
    # Box & Median
    for i in range(len(bplot['boxes'])):
        box = bplot['boxes'][i]
        median = bplot['medians'][i]
        color = colors[i]
        
        art.setp(box, color=color, linewidth=2)
        art.setp(median, color=color, linewidth=2)
        
        boxX = []
        boxY = []
        for j in box.get_xdata():
            boxX.append(j)
        for j in box.get_ydata():
            boxY.append(j)
        boxCoords = zip(boxX,boxY)
        boxPolygon = Polygon(boxCoords, facecolor=color, alpha=0.66)
        ax.add_patch(boxPolygon)
    
def plotMapBars(lOrg, title, fname, svg=False):
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
                color=['#3366CC', 'orange', '#D32626'])
        if index == 0:
            patch[0].set_label('Size')
            patch[1].set_label('Mapped to Kegg')
            patch[2].set_label('Kegg reactions')
    
    plt.xticks([0.2 + 1 * x for x in range(len(lOrg))] , [x[0] for x in lOrg])
    plt.ylim(0, maxprots + maxprots * 0.33)
    plt.title(title)
    plt.legend(loc='best',prop={'size':6})
    
    if svg:
        fname += '.svg'
    else:
        fname += '.png'
    
    plt.savefig(fname)

    logger.info('%s graph saved (%s)'%(title, fname))
    
def plotPanGenome(core, acc, uni, svg=False):
    plt.clf()
    colors=('#D32626','#3366CC','#33CC33')
    patches = plt.pie([core, acc, uni], colors=colors,
                      explode=(0.1,0.01,0.01),
                      autopct='%1.1f%%',
                      shadow=True)
    plt.legend((patches[0][0], patches[0][1], patches[0][2]),
               ('Core','Accessory','Unique'),
               loc=(0,-.1),prop={'size':6})
    plt.title('PanGenome shape')
    
    if svg:
        fname = 'pangenome_shape.svg'
    else:
        fname = 'pangenome_shape.png'
    
    plt.savefig(fname)
    
    logger.info('PanGenome shape graph saved (%s)'%fname)

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
