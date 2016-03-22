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
import math
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
        if not create.create():
            return False
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
    
    # Add to the project
    biolog = Biolog(project)
    
    if not os.path.exists(filename):
        logger.error('Phenomic file %s may not be present'%(filename))
        return False
    
    org = Organism(project)
    if not org.isOrg(orgID):
        logger.warning('Organism %s is not present yet!'%orgID)
        return False
    
    filename = os.path.abspath(filename)
    
    bparser = BiologParser(filename, 
                           set([x.plate_id
                                for x in biolog.getPlates()]))
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
    orgFound = True
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
        orgFound = False
    
    # Prepare a series of Plate objects to catch the replicas
    # (replicas will be handled by the db tough)
    dPlates={}
    for plate in bparser.plates:
        # Check if some plateIDs are unknown
        if not biolog.isPlate(plate.plate_id):
            logger.warning('Plate ID (%s) not present in the project, skipping this plate'%plate.plate_id)
            logger.warning('Or you can import your custom plate with the import-plates command')
            bparser.plates.remove(plate)
            continue
        #
        if plate.plate_id not in dPlates:
            dPlates[plate.plate_id] = Plate(plate.plate_id)
        if plate.strain in dPlates[plate.plate_id].strains:
            plate.replica = len(dPlates[plate.plate_id].strains[plate.strain]) + 1
        else:
            plate.replica = 1
        dPlates[plate.plate_id].addData(plate.strain, plate)
    
    # Grep the wells
    wells = [w for plate in dPlates.itervalues() for w in plate.getWells()]
    
    # Manually add the OrgID
    if not orgFound:
        for w in wells:
            w.strain = orgID
    
    biolog.addWells(wells, clustered=False)
    # If we have parsed a yaml/json we may have the parameters as well
    biolog.addWells(wells, clustered=True, imported=True)
    
    logger.info('Added phenome %s, having %d biolog plates (%d wells)'%
                (orgID, len(bparser.plates), len(wells)))
    
    return True

def dPhenomeMultiAdd(project, filename):
    '''
    Add a single phenomic file with multiple organisms in it
    '''
    from ductape.phenome.biolog import BiologParser, Plate
    
    if not os.path.exists(filename):
        logger.error('Phenomic file %s may not be present'%(filename))
        return False
    
    # Add to the project
    biolog = Biolog(project)
    
    filename = os.path.abspath(filename)
    
    bparser = BiologParser(filename, set([x.plate_id
                                          for x in biolog.getPlates()]))
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
                # Check if some plateIDs are unknown
                if not biolog.isPlate(plate.plate_id):
                    logger.warning('Plate ID (%s) not present in the project, skipping this plate'%plate.plate_id)
                    logger.warning('Or you can import your custom plate with the import-plates command')
                    bparser.plates.remove(plate)
                    continue
                #
                if plate.plate_id not in dPlates:
                    dPlates[plate.plate_id] = Plate(plate.plate_id)
                dPlates[plate.plate_id].addData(plate.strain, plate)
                if plate.strain in dPlates[plate.plate_id].strains:
                    plate.replica = len(dPlates[plate.plate_id].strains[plate.strain]) + 1
                else:
                    plate.replica = 1
        
        # Grep the wells
        wells = [w for plate in dPlates.itervalues() 
                 for w in plate.getWells()]
        
        biolog.addWells(wells, clustered=False)
        # If we have parsed a yaml/json we may have the parameters as well
        biolog.addWells(wells, clustered=True, imported=True)
        
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
    from ductape.phenome.biolog import BiologZero
    
    biolog = Biolog(project)
    
    if not isPhenome(project):
        logger.warning('No phenotypic data available!')
        return False
    
    if blankfile:
        logger.info('Going to use a blank file for zero subtraction')
        
        if not os.path.exists(blankfile):
            logger.error('Blank file %s may not be present'%(blankfile))
            return False
        
        bparser = BiologParser(blankfile, set([x.plate_id
                                               for x in biolog.getPlates()]))
        bparser.parse()
        
        if len(bparser.plates) == 0:
            logger.warning('The blank file contains no plates!')
            return False
    
    zeroPlates = [x.plate_id for x in biolog.getZeroSubtractablePlates()]
    
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
        # We have to provide the plates which can be subtracted,
        # along with the information on the control well of each well
        controlWells = {}
        for plate_id, zero_well_id in biolog.getControlWells():
            controlWells[plate_id] = controlWells.get(plate_id, set())
            controlWells[plate_id].add(zero_well_id)
        
        zeroWells = {}
        for plate_id, well_id, zero_well_id in biolog.getControlPairs():
            zeroWells[plate_id] = zeroWells.get(plate_id, {})
            zeroWells[plate_id][well_id] = zero_well_id
        
        zsub = BiologZero(plates, zeroPlates=zeroPlates,
                          controlWells=controlWells,
                          zeroWells=zeroWells)
        
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
        logger.warning('The parameters and the activity must be recalculated')
        biolog.delWellsParams(wells)
    
    return True

def dPhenomeTrim(project, trimTime = None):
    '''
    Takes all the biolog data available and performs the signals trimming
    The minimum time available will be used as global time maximum for all signals
    To be used when there are very different ending times inside an experiment

    If trimTime is provided, that time will be used for the trim
    '''
    from ductape.phenome.biolog import getPlates, Experiment
    
    biolog = Biolog(project)
    
    if not isPhenome(project):
        logger.warning('No phenotypic data available!')
        return False
    
    sigs = [s for s in biolog.getAllSignals()]
    plates = [p for p in getPlates(sigs, nonmean=True)]
    isZero = biolog.atLeastOneZeroSubtracted()

    if len(plates) == 0:
        logger.warning('No phenomic data available, skipping trimming')
        return True
    else:
        logger.info('Trimming %d phenomic plates'%len(plates))
    
    if trimTime is not None:
        logger.info('Trimming plates at %f'%trimTime)

    zeroPlates = [x.plate_id for x in biolog.getZeroSubtractablePlates()]
    
    exp = Experiment(plates=plates, zero=isZero, zeroPlates=zeroPlates)
    mtime = exp.trim(trimTime)
    
    logger.info('Trimmed %d plates at %f'%(len(plates), mtime))
    
    logger.info('Updating the plates')
    # Add to the project
    biolog = Biolog(project)
    biolog.updateSignals(exp.getWells(False))
    
    logger.warning('The parameters and the activity must be recalculated')
    biolog.delWellsParams(exp.getWells(False))
        
    return True

def dPhenomePurge(project, policy, delta=1, filterplates=[], replica=None):
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

    zeroPlates = [x.plate_id for x in biolog.getZeroSubtractablePlates()]
    
    exp = Experiment(plates=plates, zero=isZero, zeroPlates=zeroPlates)
    
    if delta >= exp.getMaxActivity():
        logger.warning('The delta activity threshold is higher than the maximum '+
                       'activity found (%d vs. %d)'%(delta, exp.getMaxActivity()))
    
    if not exp.purgeReplicas(delta=delta,policy=policy,replica=replica):
        logger.error('Could not purge the phenomic experiment')
        return False

    # Move the discarded wells
    biolog.moveDiscardedWells(exp.discarded)
    logger.info('Purged %d phenomic experiments'%len(exp.discarded))
    
    return True

def dPhenomeRestore(project, plates=[], replica=None):
    biolog = Biolog(project)
    
    if not biolog.atLeastOnePurged():
        logger.warning('No phenomic experiment to be restored')
        return True
    
    howmany = biolog.restoreDiscardedWells(plates, replica)
    
    logger.info('Restored %d phenomic experiments'%howmany)
    
    return True

def dGenomeAnnotate(project, noWrite=False):
    # Which project are we talking about?
    kind = dSetKind(project)
    
    proj = Project(project)
    proj.getProject()
    
    if kind != 'pangenome':
        logger.warning("Kegg annotation merge can be performed only with pangenome analysis")
        return False
    
    if not proj.isPanGenome():
        logger.warning("The pangenome has not been predicted yet")
        return False
    
    # Start to merge the annotations
    merged, mergedg, multiple = mergeKegg(project)
    
    if not noWrite:
        genome = Genome(project)
        genome.addKOs(merged, True)
    
        logger.info('Merged %d KEGG annotations'%len(merged))
        logger.info('%d orthologous groups have been re-annotated'%len(mergedg))
        
    else:
        logger.info('%d KEGG annotations can be added'%len(merged))
        logger.info('%d orthologous groups may be re-annotated'%len(mergedg))
        
    if len(multiple) >= 1:
        logger.warning('Found %d orthologous groups with more than one KO entry'%len(multiple))
                     
    return True

def dGenomeDeAnnotate(project):
    howmany = Genome(project).howManyMergedKOs()

    if howmany == 0:
        logger.info('No merged KO links have been found')
        return True
    
    Genome(project).delMergedKOs()
    
    now = Genome(project).howManyMergedKOs()
    if now > 0:
        logger.error('%d merged KO links are still present!'%now)
        return False

    logger.info('%d merged KO links have been removed'%howmany)
    logger.warning('You may want to re-run some analysis')
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
        f = open('single_stats.tsv', 'w')
        
        # Header
        header = '\t'.join( ['ID', 'name', 'description', 'proteome size',
                                'mapped to kegg', 'KEGG orthology IDs',
                                'pathways', 'reactions', 'unique reactions',
                                'exclusive reaction IDs'] )
        if doPrint:
            print(header)
        else:
            logger.info(header)
        f.write(header+'\n')
        
        eReacts = kegg.getExclusiveReactions()
        
        lOrg = []
        for org in organism.getAll():
            org_id = org.org_id
            name = org.name if org.name else 'NONE'
            description = org.description if org.description else 'NONE'
            
            prots = genome.howMany(org_id)
            
            mapped, ko, react, path, unireact, ereact = (kegg.howManyMapped(org_id),
                                        kegg.howManyKO(org_id),
                                        kegg.howManyReactions(org_id),
                                        kegg.howManyPathways(org_id),
                                        kegg.howManyUniqueReactions(org_id),
                                        len(eReacts[org_id]))
            
            stats = '\t'.join( [str(x) for x in [org_id, name, description,
                                                 prots, mapped, ko, path,
                                                 react, unireact, ereact]] )
            if doPrint:
                print(stats)
            else:
                logger.info(stats)
            f.write(stats+'\n')
                
            lOrg.append([org_id, prots, mapped, react, unireact])
            
        f.close()
        logger.info('Table also saved in file %s'%('single_stats.tsv'))       
        plotMapBars(lOrg, 'Single genomes statistics', 'single', svg)
        
        if proj.isPanGenome():
            logger.info('Pangenome stats (orthologs)')
            logger.warning('Please note that the dispensable genome includes'+
                           ' the accessory and unique genome')
            
            # Pangenome stats
            f = open('pangenome_stats.tsv', 'w')
            
            # Header
            header = '\t'.join( ['kind', 'size',
                                'mapped to kegg', 'KEGG orthology IDs',
                                'pathways', 'reactions', 'unique reactions', 
                                'exclusive reaction IDs'] )
            if doPrint:
                print(header)
            else:
                logger.info(header)
            f.write(header+'\n')
                
            core, disp, acc, uni = (genome.getLenCore(), genome.getLenDisp(),
                              genome.getLenAcc(), genome.getLenUni())

            ecore, edisp, eacc, euni = kegg.getExclusiveReactionsPanGenome()

            stats = []
            stats.append('\t'.join( [str(x) for x in ['core', core,
                                 kegg.howManyMapped(pangenome='core'),
                                 kegg.howManyKO(pangenome='core'),
                                 kegg.howManyPathways(pangenome='core'),
                                 kegg.howManyReactions(pangenome='core'),
                                 kegg.howManyUniqueReactions(pangenome='core'),
                                 len(ecore)]]))
            stats.append('\t'.join( [str(x) for x in ['dispensable', disp,
                                 kegg.howManyMapped(pangenome='dispensable'),
                                 kegg.howManyKO(pangenome='dispensable'),
                                 kegg.howManyPathways(pangenome='dispensable'),
                                 kegg.howManyReactions(pangenome='dispensable'),
                                 kegg.howManyUniqueReactions(pangenome='dispensable'),
                                 len(edisp)]]))
            stats.append('\t'.join( [str(x) for x in ['accessory', acc,
                                 kegg.howManyMapped(pangenome='accessory'),
                                 kegg.howManyKO(pangenome='accessory'),
                                 kegg.howManyPathways(pangenome='accessory'),
                                 kegg.howManyReactions(pangenome='accessory'),
                                 kegg.howManyUniqueReactions(pangenome='accessory'),
                                 len(eacc)]]))
            stats.append('\t'.join( [str(x) for x in ['unique', uni,
                                 kegg.howManyMapped(pangenome='unique'),
                                 kegg.howManyKO(pangenome='unique'),
                                 kegg.howManyPathways(pangenome='unique'),
                                 kegg.howManyReactions(pangenome='unique'),
                                 kegg.howManyUniqueReactions(pangenome='unique'),
                                 len(euni)]]))
            
            for stat in stats:
                if doPrint:
                    print(stat)
                else:
                    logger.info(stat)
                f.write(stat+'\n')
                
            f.close()
            logger.info('Table also saved in file %s'%('pangenome_stats.tsv'))
            
            lPanGenome = [['Core', core, kegg.howManyMapped(pangenome='core'),
                           kegg.howManyReactions(pangenome='core'),
                           len(ecore)],
                          ['Dispensable', disp,
                           kegg.howManyMapped(pangenome='dispensable'),
                           kegg.howManyReactions(pangenome='dispensable'),
                           len(edisp)],
                          ['Accessory', acc,
                           kegg.howManyMapped(pangenome='accessory'),
                           kegg.howManyReactions(pangenome='accessory'),
                           len(eacc)],
                          ['Unique', uni,
                           kegg.howManyMapped(pangenome='unique'),
                           kegg.howManyReactions(pangenome='unique'),
                           len(euni)]]
 
            plotMapBars(lPanGenome, 'PanGenome statistics', 'pangenome_stats',
                        svg, labels=['Size', 'Mapped to Kegg',
                                'Kegg reactions', 'Excusive Kegg reaction IDs'])
            plotPanGenome(core, acc, uni, svg)
            
            logger.info('Pangenome stats (reactions)')
            logger.warning('Please note that here we consider the presence of '+
                           'distinct reaction IDs in each organism')
            # Pangenome stats
            f = open('pangenome_reactions_stats.tsv', 'w')
            
            # Header
            header = '\t'.join( ['kind', 'distinct reaction IDs'] )
            if doPrint:
                print(header)
            else:
                logger.info(header)
            f.write(header+'\n')
                
            conserved = set([r.re_id for r in kegg.getConservedReactions()])
            variable = set([r.re_id for r in kegg.getVariableReactions()]) 

            stats = []
            stats.append('\t'.join( [str(x) for x in ['conserved',
                                                      len(conserved)]]))
            stats.append('\t'.join( [str(x) for x in ['variable',
                                                      len(variable)]]))
            
            for stat in stats:
                if doPrint:
                    print(stat)
                else:
                    logger.info(stat)
                f.write(stat+'\n')
                    
            f.close()
            logger.info('Table also saved in file %s'%(
                                               'pangenome_reactions_stats.tsv'))
            plotPanGenomeReactions(len(conserved), len(variable), svg)
    
    elif kind == 'mutants':
        refs = [org.org_id
                    for org in organism.getAll()
                    if not organism.isMutant(org.org_id)]
        
        # Header
        header = '\t'.join( ['ID', 'name', 'description', 'kind', 'proteome size',
                                'mapped to kegg', 'reactions', 'unique reactions',
                                'exclusive reaction IDs'] )
        
        logger.warning('Please note that for deletion mutants the '+
                       '"exclusive reaction IDs" column refers to the ones '+
                       'that are present in the wild-type and NOT '+
                       'in the mutant')
        for ref_id in refs:
            logger.info('Mutants of %s stats'%ref_id)
            
            f = open('%s_stats.tsv'%ref_id, 'w')
            
            if doPrint:
                print(header)
            else:
                logger.info(header)
            f.write(header+'\n')
            
            muts = [x for x in organism.getOrgMutants(ref_id)]
            
            eReacts = kegg.getExclusiveReactionsMutants(ref_id, muts)
            
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
                
                mapped, react, unireact, ereact = (kegg.howManyMapped(org_id),
                                kegg.howManyReactions(org_id),
                                kegg.howManyUniqueReactions(org_id),
                                len(eReacts[org_id]))
        
                if mkind == 'deletion':
                    mapped = kegg.howManyMapped(ref_id) - mapped
                    react = kegg.howManyReactions(ref_id) - react
                    # TODO: not sure this calculation is correct
                    unireact = kegg.howManyUniqueReactions(ref_id) - unireact
                    #
                    
                elif mkind == 'insertion':
                    mapped += kegg.howManyMapped(ref_id)
                    react += kegg.howManyReactions(ref_id)
                    # TODO: not sure this calculation is correct
                    unireact += kegg.howManyUniqueReactions(ref_id)
                    #
                
                stats = '\t'.join( [str(x) for x in [org_id, name, description,
                                                 mkind, prots, mapped,
                                                 react, unireact, ereact]] )
                if doPrint:
                    print(stats)
                else:
                    logger.info(stats)
                f.write(stats+'\n')
                
                lOrg.append([org_id, prots, mapped, react, ereact])
        
            f.close()
            logger.info('Table also saved in file %s'%(
                                               '%s_stats.tsv'%ref_id))
            
            plotMapBars(lOrg, 'Wild-type (%s) and mutants statistics'%ref_id,
                        '%s'%ref_id, svg,
                        labels=['Size', 'Mapped to Kegg',
                                'Kegg reactions', 'Excusive Kegg reaction IDs'])
    
    else:
        logger.info('No statistics can be computed at this time')
        return False

    return True

def dPhenomeStats(project, activity=5, delta=3, svg=False, doPrint=True):
    from ductape.phenome.biolog import getPlates, Experiment
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
    
    zeroPlates = [x.plate_id for x in biolog.getZeroSubtractablePlates()]
    
    exp = Experiment(plates=plates, zero=isZero,
                     category=category, categorder=categorder,
                     zeroPlates=zeroPlates)
    
    exp.plot(svg=svg)
    
    # Max value for activity
    maxAct = exp.getMaxActivity()
    
    # Check the activity threshold
    if activity > maxAct:
        logger.warning('The activity threshold is higher than the maximum '+
                       'activity found (%d vs. %d)'%(activity, maxAct))
        return False
        
    # Check the activity delta value
    if delta >= maxAct:
        logger.warning('The delta activity threshold is higher than the maximum '+
                       'activity found (%d vs. %d)'%(delta, maxAct))
        return False
    
    ############################################################################
    # Activity distribution
    
    logger.info('Activity distributions')
    
    # Fake plot top get the correct bin centers
    y,binEdges=np.histogram(range(10),bins=maxAct + 1)
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
            y,binEdges=np.histogram(x,bins=maxAct + 1)
            ax.plot(bincenters,y,'-o', color='black', linewidth=2)
            
            ax.set_ylim(0,max([x for x in d.itervalues()]) + 20)
            ax.set_xlim(min(bincenters)-1,max(bincenters)+1)
            
            x0,x1 = ax.get_xlim()
            y0,y1 = ax.get_ylim()
            ax.set_aspect((x1-x0)/(y1-y0))
            
            ax.set_xticks(bincenters)
            ax.set_xticklabels([str(x) for x in range(maxAct + 1)])
            
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
                
        y,binEdges=np.histogram(x,bins=maxAct + 1)
        ax.plot(bincenters,y,'-o', color=dcolors[org.org_id], linewidth=2,
                alpha=0.66, label=org.org_id)
        
        maxv.append(max([x for x in d.itervalues()]))
        
    ax.set_ylim(0,max(maxv) + 20)
    ax.set_xlim(min(bincenters)-1,max(bincenters)+1)
    
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    ax.set_aspect((x1-x0)/(y1-y0))
    
    ax.set_xticks(bincenters)
    ax.set_xticklabels([str(x) for x in range(maxAct + 1)])
    
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
            
            y,binEdges=np.histogram(x,bins=maxAct + 1)
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
        ax.set_xticklabels([str(x) for x in range(maxAct + 1)])
            
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
            
            y,binEdges=np.histogram(x,bins=maxAct + 1)
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
        ax.set_xticklabels([str(x) for x in range(maxAct + 1)])
            
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
            
            ax.set_xlim(-1,maxAct + 1)
            ax.set_xticks(range(maxAct + 1))
            
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
        ax.set_xlim(-1,maxAct + 1)
        ax.set_xticks(range(maxAct + 1))
        
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
            ax.set_xlim(-1,maxAct + 1)
            ax.set_xticks(range(maxAct + 1))
            
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
    
    f = open('active_stats.tsv', 'w')
    
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
        print(header)
    else:
        logger.info(header)
    f.write(header+'\n')
    
    header = '\t'.join( ['Category'] + orgs )
    if doPrint:
        print(header)
    else:
        logger.info(header)
    f.write(header+'\n')
        
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
            print(line)
        else:
            logger.info(line)
        f.write(line+'\n')
    
    f.close()
    logger.info('Table also saved in file %s'%('active_stats.tsv'))
    
    if kind == 'single':
        return True
    
    logger.info('Active differences stats')
    f = open('active_diffs_stats.tsv', 'w')
    
    header = '\t'.join( ['Category', 'Average difference',
             'Main differences (%% of wells whose average difference >= %d)'%delta] )
    if doPrint:
        print(header)
    else:
        logger.info(header)
    f.write(header+'\n')
        
    for categ in categorder:
        pwlist = filter(lambda x: x[0] in category[categ],
                    [pw for pw in getOrder(project)])
        
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
            print(line)
        else:
            logger.info(line)
        f.write(line+'\n')
    
    f.close()
    logger.info('Table also saved in file %s'%('active_diffs_stats.tsv'))
       
    ############################################################################
    
    if kind == 'single':
        return True
    
    logger.info('Unique metabolic functions')
    f = open('unique_stats.tsv', 'w')
    
    header = 'Unique metabolic functions (%% of wells with delta >= %d)'%delta
    if doPrint:
        print(header)
    else:
        logger.info(header)
    f.write(header+'\n')
    
    header = '\t'.join( [''] + orgs )
    if doPrint:
        print(header)
    else:
        logger.info(header)
    f.write(header+'\n')
    
    UpUnique = {}
    for oid in orgs:
        UpUnique[oid] = 0
    DownUnique = {}
    for oid in orgs:
        DownUnique[oid] = 0
    pwcount = 0
    for pid, wid in getOrder(project):
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
        print(line)
    else:
        logger.info(line)
    f.write(line+'\n')
    
    line = '\t'.join(['Less active'] + [str(DownUnique[oid]) for oid in orgs])    
    if doPrint:
        print(line)
    else:
        logger.info(line)
    f.write(line+'\n')
            
    line = '\t'.join(['Total'] + [str(UpUnique[oid] + DownUnique[oid]) for oid in orgs])    
    if doPrint:
        print(line)
    else:
        logger.info(line)
    f.write(line+'\n')
    
    line = '\t'.join(['%'] + [str(((UpUnique[oid] + DownUnique[oid])/total)*100) for oid in orgs])    
    if doPrint:
        print(line)
    else:
        logger.info(line)
    f.write(line+'\n')
    
    f.close()
    logger.info('Table also saved in file %s'%('unique_stats.tsv'))
    
    return True

def dPhenomeRings(project, delta=1, difforg=None, svg=False,
        param='activity'):
    from ductape.phenome.biolog import getPlates, Experiment
    
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
    
    zeroPlates = [x.plate_id for x in biolog.getZeroSubtractablePlates()]
    
    exp = Experiment(plates=plates, zero=isZero,
                     category=category, categorder=categorder,
                     zeroPlates=zeroPlates)
    
    # Max value for the chosen parameter
    maxAct = exp.getMaxParam(param)
    
    # Check the activity delta value
    # We don't check the delta for other parameters
    # Hopefully they have an higher range
    if param == 'activity' and delta >= maxAct:
        logger.warning('The delta activity threshold is higher than the maximum '+
                       'activity found (%d vs. %d)'%(delta, maxAct))
        return False
    
    ############################################################################
    # Activity rings (!!!)
    # Thanks to stackoverflow for that
    # http://stackoverflow.com/questions/12803883
    
    logger.info('Activity ring')
    if param != 'activity':
        logger.info('Using %s instead of activity'%param)
    
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
    if len(orgs) > 10:
        i = 0.05
        i_space = 0.05
        text_incr = 0.015
        i_incr = 0.07
        t_size = 10
    else:
        i = 0.1
        i_space = 0.1
        text_incr = 0.03
        i_incr = 0.15
        t_size = 17
    for org_id in orgs:
        radius = np.linspace(i, i+i_space, 10)
        theta = np.linspace(0, 2*np.pi, 628)
        R,T  = np.meshgrid(radius,theta)
    
        ax.pcolor(T, R, np.array([[1 for y in range(10)] for x in range(628)]),
                  cmap=cm.Greys,
                  vmin=0, vmax=maxAct)
        
        ax.text(0, i+text_incr, org_id, size=t_size, weight='black', alpha=0.77, ha='center')
        
        i += i_incr
        
    # Category and plate/well fail-safe tweaks
    # What if a plate/well is missing?
    # What is the category order?
    categpworder = {}
    for categ in categorder:
        categpworder[categ] = []
        for pid, wid in getOrder(project, sorted(category[categ])):
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
                except Exception as e:
                    logger.warning(e)
                    acts.append(np.nan)
                    continue
                if (kind == 'mutants' or difforg) and org_id in muts:
                    # Check if the reference has a value for this
                    try:
                        ref_id = muts[org_id]
                        refact = getattr(exp.sumexp[w.plate_id][w.well_id][ref_id], param)
                        if param == 'activity' and abs(getattr(w, param) - refact) <= delta:
                            acts.append(np.nan)
                        else:
                            acts.append(getattr(w, param) - refact)
                    except Exception as e:
                        logger.warning(e)
                        acts.append(np.nan)
                else:
                    if getattr(w, param) is None:
                        acts.append(np.nan)
                    else:
                        acts.append(getattr(w, param))
                    
        radius = np.linspace(i, i+0.2, 10)
        theta = np.linspace(0, 2*np.pi, len(acts))
        R,T  = np.meshgrid(radius,theta)
        
        if (kind == 'mutants' or difforg) and org_id in muts:
            cmap = cm.PuOr
            vmin = -maxAct
            vmax = maxAct
        else:
            cmap = cm.RdYlGn
            vmin = 0
            vmax = maxAct
            
        cmap.set_under('#F8F8F8',1.)
        
        ax.pcolor(T, R, np.array([[x for y in range(10)] for x in acts]),
                  cmap=cmap,
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
    
    if len(orgs) > 10:
        l_width = 15
    else:
        l_width = 35    
    for categ in categorder:
        if len(categpworder[categ]) == 0:continue
        
        # This category proportion over the others
        categprop = float(len(categpworder[categ])) / total
        
        stop_arch += 2*np.pi*categprop
        theta = np.linspace(start_arch, stop_arch, 100)
        
        ax.plot(theta, [i for t in theta], color=categcolor[categ],
                linewidth=l_width, label=categ)
        
        start_arch += 2*np.pi*categprop

    #Turn off polar labels
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    ax.set_rmax(i)

    # Colorbar
    import matplotlib.colors as colors
    cNorm  = colors.Normalize(vmin=0, vmax=maxAct)
    scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cm.RdYlGn)
    scalarMap.set_array(np.array(range(int(maxAct) + 1)))
    cax = fig.add_axes([0.93, 0.2, 0.03, 0.6])
    cax.text(0.50, 1.01, param, size=20, ha='center')
    plt.colorbar(scalarMap, cax=cax)
    
    if (kind == 'mutants' or difforg):
        cNorm  = colors.Normalize(vmin=-maxAct, vmax=maxAct)
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cm.PuOr)
        scalarMap.set_array(np.array(range(-int(maxAct),int(maxAct),20)))
        cax = fig.add_axes([0.04, 0.2, 0.03, 0.6])
        cax.text(0.50, 1.01, 'Delta %s'%param, size=20, ha='center')
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
    pkind = dSetKind(project)
    
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
        
    proj = Project(project)
    
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
    
    if proj.isPanGenome():
        merged, mergedg, multiple = mergeKegg(project)
        if len(merged) > 0:
            nprots = set([x[0] for x in merged])
            logger.warning('Annotation on %d proteins can be improved'%len(nprots))
            logger.warning('%d orthologous groups would be affected'%len(mergedg))
        else:
            for org in organism.getAll():
                fname = 'ko_merged_%s.tsv'%org.org_id
                fout = open(fname,'w')
                fout.write('#%s\t%s\n'%('prot_id', 'ko_id'))
                i = 0
                # Get the merged annotations (if any)
                for prot_id, ko_id in kegg.getAllKO(org.org_id, merged=True):
                    fout.write('%s\t%s\n'%(prot_id, ko_id.lstrip('ko:')))
                    i += 1
                fout.close()
                
                if i == 0:
                    os.remove(fname)
                    logger.info('No merged KO links available for %s'
                                        %org.org_id)
                else:
                    logger.warning('Saved %d merged KO links for %s (%s)'
                                %(i, org.org_id, fname))
        
        if len(multiple) > 0:
            logger.warning('Found %d orthologous groups with more than one KO entry'%len(multiple))
            
            fname = 'ko_multiple.tsv'
            fout = open(fname,'w')
            fout.write('#%s\n'%('group_id'))
            i = 0
            for group_id in multiple:
                fout.write('%s\n'%group_id)
                i += 1
            fout.close()
            
            logger.warning('Saved %d orthologs with multiple KO links (%s)'
                            %(i, fname))
            
            for org in organism.getAll():
                fname = 'ko_multiple_%s.tsv'%org.org_id
                fout = open(fname,'w')
                fout.write('#%s\t%s\n'%('prot_id', 'ko_id'))
                i = 0
                # Get the multiple annotations (if any)
                for prot_id, ko_id in kegg.getMultipleKOs(org.org_id):
                    fout.write('%s\t%s\n'%(prot_id, ko_id.lstrip('ko:')))
                    i += 1
                fout.close()
                
                if i == 0:
                    os.remove(fname)
                    logger.info('No mutiple KO links available for %s'
                                        %org.org_id)
                else:
                    logger.warning('Saved %d multiple KO links for %s (%s)'
                                %(i, org.org_id, fname))
        
    logger.info('Exporting Kegg reactions data')
    
    for org in organism.getAll():
        fname = 'reactions_%s.tsv'%org.org_id
        fout = open(fname,'w')
        fout.write('#%s\t%s\t%s\t%s\t%s\n'%('prot_id', 're_id', 'name', 'description', 'pathway(s)'))
        i = 0
        for re in kegg.getAllReactions(org.org_id):
            paths = ','.join([p.path_id.lstrip('path:') for p in kegg.getReactPath(re.re_id)])
            fout.write('%s\t%s\t%s\t%s\t%s\n'%(re.prot_id, re.re_id.lstrip('rn:'), re.name,
                                  re.description, paths))
            i += 1
        fout.close()
        
        if i == 0:
            os.remove(fname)
            logger.warning('No Kegg reactions available for %s'%org.org_id)
        else:
            logger.info('Saved %d Kegg reactions links for %s (%s)'%
                    (i, org.org_id, fname))
    
    if pkind != 'mutants':
        eReacts = kegg.getExclusiveReactions()
        for org in organism.getAll():
            if len(eReacts[org.org_id]) == 0:
                logger.warning('No exclusive Kegg reactions available for %s'%org.org_id)
                continue
            
            fname = 'reactions_exclusive_%s.tsv'%org.org_id
            fout = open(fname,'w')
            fout.write('#%s\t%s\t%s\t%s\n'%('re_id', 'name', 'description', 'pathway(s)'))
            i = 0
            for re_id in eReacts[org.org_id]:
                re = kegg.getReaction(re_id)
                paths = ','.join([p.path_id.lstrip('path:') for p in kegg.getReactPath(re.re_id)])
                fout.write('%s\t%s\t%s\t%s\n'%(re_id.lstrip('rn:'), re.name,
                                  re.description, paths))
                i += 1
            fout.close()
            
            logger.info('Saved %d exclusive Kegg reactions for %s (%s)'%
                        (i, org.org_id, fname))
    else:
        logger.warning('Please note that for deletion mutants the '+
                       '"exclusive reaction IDs" column refers to the ones '+
                       'that are present in the wild-type and NOT '+
                       'in the mutant')
        
        refs = [org.org_id
                    for org in organism.getAll()
                    if not organism.isMutant(org.org_id)]
        
        for ref_id in refs:
            muts = [x for x in organism.getOrgMutants(ref_id)]
            
            eReacts = kegg.getExclusiveReactionsMutants(ref_id, muts)
            
            for org_id in [ref_id] + muts:
                if len(eReacts[org_id]) == 0:
                    logger.warning('No exclusive Kegg reactions available for %s'%org_id)
                    continue
                
                fname = 'reactions_exclusive_%s.tsv'%org_id
                fout = open(fname,'w')
                fout.write('#%s\t%s\t%s\t%s\n'%('re_id', 'name', 'description', 'pathway(s)'))
                i = 0
                for re_id in eReacts[org_id]:
                    re = kegg.getReaction(re_id)
                    paths = ','.join([p.path_id.lstrip('path:') for p in kegg.getReactPath(re.re_id)])
                    fout.write('%s\t%s\t%s\t%s\n'%(re_id.lstrip('rn:'), re.name,
                                  re.description, paths))
                    i += 1
                fout.close()
                
                logger.info('Saved %d exclusive Kegg reactions for %s (%s)'%
                            (i, org_id, fname))
      
    if proj.isPanGenome():
        ecore, edisp, eacc, euni = kegg.getExclusiveReactionsPanGenome()
        preact = [('core', ecore), ('dispensable', edisp),
                  ('accessory', eacc), ('unique', euni)]
        
        for label, er in preact:
            if len(er) == 0:
                logger.warning('No exclusive Kegg reactions available'+
                               ' for %s genome'%label)
                continue
            
            fname = 'reactions_exclusive_%s.tsv'%label
            fout = open(fname,'w')
            fout.write('#%s\t%s\t%s\t%s\n'%('re_id', 'name', 'description', 'pathway(s)'))
            i = 0
            for re_id in er:
                re = kegg.getReaction(re_id)
                paths = ','.join([p.path_id.lstrip('path:') for p in kegg.getReactPath(re.re_id)])
                fout.write('%s\t%s\t%s\t%s\n'%(re_id.lstrip('rn:'), re.name,
                              re.description, paths))
                i += 1
            fout.close()
            
            logger.info('Saved %d exclusive Kegg reactions for %s genome (%s)'%
                        (i, label, fname))
            if label == 'dispensable':
                logger.warning('Please note that the dispensable genome includes'+
                           ' the accessory and unique genome')
        
        ecore = [r for r in kegg.getConservedReactions()]
        edisp = [r for r in kegg.getVariableReactions()]
        preact = [('conserved', ecore), ('variable', edisp)]
        
        for label, er in preact:
            if len(er) == 0:
                logger.warning('No Kegg reactions available'+
                               ' for %s reactions set'%label)
                continue
            fname = 'reactions_%s.tsv'%label
            fout = open(fname,'w')
            fout.write('#%s\t%s\t%s\t%s\n'%('re_id', 'name', 'description', 'pathway(s)'))
            i = 0
            for r in er:
                re_id = r.re_id
                re = kegg.getReaction(re_id)
                paths = ','.join([p.path_id.lstrip('path:') for p in kegg.getReactPath(re.re_id)])
                fout.write('%s\t%s\t%s\t%s\n'%(re_id.lstrip('rn:'), re.name,
                              re.description, paths))
                i += 1
            fout.close()
            
            logger.info('Saved %d exclusive Kegg reactions for %s genome (%s)'%
                        (i, label, fname))
        
    logger.info('Exporting EC numbers')
    
    for org in organism.getAll():
        fname = 'ecnumbers_%s.tsv'%org.org_id
        fout = open(fname,'w')
        fout.write('#%s\t%s\n'%('prot_id', 'EC_number'))
        i = 0
        for prot_id, ec_ids in kegg.getAllECNumbers(org.org_id):
            while '  ' in ec_ids:
                ec_ids = ec_ids.replace('  ', ' ')
            for ec_id in ec_ids.split():
                fout.write('%s\t%s\n'%(prot_id, ec_id))
                i += 1
        fout.close()
        
        if i == 0:
            os.remove(fname)
            logger.warning('No EC numbers available for %s'%org.org_id)
        else:
            logger.info('Saved %d EC numbers links for %s (%s)'%
                    (i, org.org_id, fname))
    
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
                                           '\t'.join(dG[group.group_id])))
            for group in genome.getAcc():
                fout.write('%s\t%s\t%s\n'%(group.group_id,
                                           'accessory',
                                           '\t'.join(dG[group.group_id])))
            for group in genome.getUni():
                fout.write('%s\t%s\t%s\n'%(group.group_id,
                                           'unique',
                                           '\t'.join(dG[group.group_id])))
            fout.close()
            
            logger.info('Exported orthologs informations (%s)'%fname)
    
    return True

def dBiologImport(project, infile):
    biolog = Biolog(project)
    
    logger.info('Importing custom plate(s)')
    
    before = 0
    for r in biolog.getPlates():
        before += 1
    
    try:
        biolog.importBiolog(infile)
    except Exception as e:
        logger.error('Could not import the custom plate(s)!')
        logger.error(e)
        return False
    
    after = 0
    for r in biolog.getPlates():
        after += 1
        
    logger.info('Imported %d custom plate(s)'%(after-before))
    
    return True

def dPhenomeExport(project, json=False):
    from ductape.phenome.biolog import getSinglePlatesFromSignals
    from ductape.phenome.biolog import toYAML, toJSON
    from ductape.common.utils import safeSubtraction
    
    biolog = Biolog(project)    
    
    logger.info('Exporting Biolog plate data')
    fname = 'biolog.tsv'
    fout = open(fname,'w')
    
    i = 0
    for row in biolog.exportBiolog():
        fout.write(row + '\n')
        i += 1
        
    fout.close()
    
    logger.info('Exported %d Biolog wells (%s)'%(i, fname))
    
    # Is there something to be exported?
    organism = Organism(project)
    
    if organism.howMany() == 0:
        logger.info('No phenomic data can be exported at this time')
        return True
    
    # Which project are we talking about?
    kind = dSetKind(project)
    
    logger.info('Exporting phenomic data for other programs')
    
    sigs = [x for x in biolog.getAllSignals()]
    for plate in getSinglePlatesFromSignals(sigs):
        logger.info('Exporting plate %s, strain %s, replica %d'%(plate.plate_id,
                                                                 plate.strain,
                                                                plate.replica))
        if not json:
            fout = open('%s_%s_%s.yml'%(plate.plate_id, plate.strain,
                                    plate.replica), 'w')
            fout.write(toYAML(plate))
            fout.close()
        else:
            fout = open('%s_%s_%s.json'%(plate.plate_id, plate.strain,
                                    plate.replica), 'w')
            fout.write(toJSON(plate))
            fout.close()
    
    logger.info('Exporting single organism(s) phenomic data')
    
    for org in organism.getAll():
        fname = 'phenome_%s.tsv'%org.org_id
        fout = open(fname,'w')
        fout.write('#' + '\t'.join(['plate_id', 'well_id', 'chemical',
                                'category',
                                'moa', 'co_id', 'replica', 'activity',
                                'min', 'max', 'height', 'plateau', 'slope',
                                'lag', 'area', 'source']) + '\n')
        i = 0
        for w in biolog.getOrgWells(org.org_id):
            wdet = biolog.getWell(w.plate_id, w.well_id)
            fout.write('\t'.join([xstr(x) for x in [w.plate_id, w.well_id, wdet.chemical,
                                  wdet.category, wdet.moa, wdet.co_id]] +
                                  [xstr(x) for x in [w.replica, w.activity,
                                                    w.min, w.max, w.height,
                                                    w.plateau, w.slope,
                                                    w.lag, w.area,
                                                    w.source]])
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
                                    'moa', 'co_id',
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
                                                [xstr(safeSubtraction(ref_act, 
                                                      biolog.getAvgActivity(
                                                                    w.plate_id,
                                                                    w.well_id,
                                                                    x)))
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

def getTotalNet(project, path_id=None):
    '''
    Get the overall Kegg metabolic net
    '''
    if path_id:
        logger.debug('Building total metabolic network (%s)'%path_id)
    else:
        logger.debug('Building total metabolic network')
    
    from ductape.kegg.net import MetabolicNet
    
    kegg = Kegg(project)

    return MetabolicNet(kegg.getAllCompounds(path_id),
                        kegg.getAllRPairsReacts(path_id=path_id))

def getOrgNet(project, org_id, path_id=None, category=None):
    '''
    Get the overall Kegg metabolic net
    '''
    if path_id:
        logger.debug('Building total metabolic network for %s (%s)'%(org_id,
                                                                     path_id))
    else:
        logger.debug('Building total metabolic network for %s'%org_id)
        
    from ductape.kegg.net import MetabolicNet, Compound
    
    kegg = Kegg(project)
    
    net = MetabolicNet(kegg.getAllCompounds(path_id),
                       kegg.getAllRPairsReacts(org_id, path_id))
    
    if category:
        logger.debug('Fetching metabolic activity for category %s'%category)
        
        biolog = Biolog(project)
        vmax = biolog.getMaxActivity()

        corg = {}
        
        # Filter by path?
        path_co = []
        if path_id:
            path_co = [x.co_id for x in kegg.getPathComps(path_id)]
            wells = [w for w in biolog.getAllCoByCateg(category)
                        if 'cpd:'+w.co_id in path_co]
        else:
            wells = [w for w in biolog.getAllCoByCateg(category)]
        for well in wells:
            act = biolog.getAvgActivity(well.plate_id, well.well_id, org_id)
            if act is not None:
                # Some co_ids are present more than once
                if well.co_id not in corg:
                    corg[well.co_id] = []
                corg[well.co_id].append(act)
    
        toremove = set()
        for k, v in corg.iteritems():
            mean = np.array(v).mean()
            if not math.isnan(float(mean)):
                corg[k] = mean
            else:
                toremove.add(k)
        for k in toremove:
            del corg[k]
        
        compounds = [Compound('cpd:'+k,kegg.getCompound('cpd:'+k).name,v,vmax) for k,v in corg.iteritems()]
        net.addNodes(compounds)
        logger.debug('Added %d metabolic activities'%len(compounds))
        
    return net

def getMutNet(project, mut_id, mut_rpairs, path_id=None, category=None):
    '''
    Get the overall Kegg metabolic net
    '''
    if path_id:
        logger.debug('Building total metabolic network for %s (%s)'%(mut_id,
                                                                     path_id))
    else:
        logger.debug('Building total metabolic network for %s'%mut_id)
        
    from ductape.kegg.net import MetabolicNet, Compound
    
    kegg = Kegg(project)
    
    net = MetabolicNet(kegg.getAllCompounds(path_id),
                       mut_rpairs)
    
    if category:
        logger.debug('Fetching metabolic activity for category %s'%category)
        
        biolog = Biolog(project)
        vmax = biolog.getMaxActivity()

        corg = {}
        
        # Filter by path?
        path_co = []
        if path_id:
            path_co = [x.co_id for x in kegg.getPathComps(path_id)]
            wells = [w for w in biolog.getAllCoByCateg(category)
                        if 'cpd:'+w.co_id in path_co]
        else:
            wells = [w for w in biolog.getAllCoByCateg(category)]
        for well in wells:
            act = biolog.getAvgActivity(well.plate_id, well.well_id, mut_id)
            if act is not None:
                # Some co_ids are present more than once
                if well.co_id not in corg:
                    corg[well.co_id] = []
                corg[well.co_id].append(act)
    
        toremove = set()
        for k, v in corg.iteritems():
            mean = np.array(v).mean()
            if not math.isnan(float(mean)):
                corg[k] = mean
            else:
                toremove.add(k)
        for k in toremove:
            del corg[k]
        
        compounds = [Compound('cpd:'+k,kegg.getCompound('cpd:'+k).name,v,vmax) for k,v in corg.iteritems()]
        net.addNodes(compounds)
        logger.debug('Added %d metabolic activities'%len(compounds))
        
    return net

def getPanGenomeNet(project, dpangenome, pangenome='all', path_id=None, category=None):
    '''
    Get a pangenomic slice of the metabolic network
    dpangenome is an input dictionary with pangenomic slices
    '''
    if pangenome not in ['all', 'conserved', 'variable']:
        logger.warning('Unknown pangenomic slice! (%s)'%pangenome)
        return None
    
    if category and pangenome != 'all':
        logger.warning('Phenomic data can be used only with all the pangenome')
        return None
    
    if path_id:
        logger.debug('Building total metabolic network for %s (%s)'%(pangenome,
                                                                     path_id))
    else:
        logger.debug('Building total metabolic network for %s'%pangenome)
        
    from ductape.kegg.net import MetabolicNet, Compound
    from itertools import combinations
    
    kegg = Kegg(project)
    
    net = MetabolicNet(kegg.getAllCompounds(path_id),
                       dpangenome[pangenome])
    
    if category:
        logger.debug('Fetching metabolic activity for category %s'%category)
        
        biolog = Biolog(project)
        vmax = biolog.getMaxActivity()
        
        corg = {}
        
        # Filter by path?
        path_co = []
        if path_id:
            path_co = [x.co_id for x in kegg.getPathComps(path_id)]
            wells = [w for w in biolog.getAllCoByCateg(category)
                        if 'cpd:'+w.co_id in path_co]
        else:
            wells = [w for w in biolog.getAllCoByCateg(category)]
        
        for well in wells:
            acts = [x.avgact 
                    for x in
                    biolog.getAvgActivityEachOrg(well.plate_id,
                                                 well.well_id)]
            if len(acts) <= 1:
                continue
            
            diffs = []
            for a, a1 in combinations(acts, 2):
                diffs.append( abs(a - a1) )
                
            avgdiff = np.array(diffs).mean()
            
            # Some co_ids are present more than once
            if well.co_id not in corg:
                corg[well.co_id] = []
            corg[well.co_id].append(avgdiff)
        
        toremove = set()
        for k, v in corg.iteritems():
            mean = np.array(v).mean()
            if not math.isnan(float(mean)):
                corg[k] = mean
            else:
                toremove.add(k)
        for k in toremove:
            del corg[k]
        
        compounds = [Compound('cpd:'+k,kegg.getCompound('cpd:'+k).name,v,vmax) for k,v in corg.iteritems()]
        net.addNodes(compounds)
        logger.debug('Added %d metabolic activities'%len(compounds))
        
    return net

def writeNet(net, path, name):
    '''
    Save a network as a gml file in the desired location
    '''
    import networkx as nx
    try:
        nx.write_gml(net.net, os.path.join(path,name),
                     nx.readwrite.gml.literal_stringizer)
    except AttributeError:
        # old version of networkx
        nx.write_gml(net.net, os.path.join(path,name))
    logger.debug('Saved network %s to %s'%(name, path))

def writeCombinedPanGenome(dvalues):
    '''
    Write down the table with the combined pangenome data
    '''
    fname = 'combined_pangenome.tsv'
    fout = open(fname,'w')
    
    fout.write('# Genomic variability Vs. Phenomic variability\n')
    fout.write('# disp/total = number of variable reaction IDs / '+
               'number of total reaction IDs \n')
    fout.write('# diffAV = mean AV difference\n')
    fout.write('# (-1, no activity available)\n')
    fout.write('\t'.join( ['pathway', 'name', 'category', 'disp/core', 'diffAV'] )
                    + '\n')
    
    # Order the pathway list for genomic variability
    allv = []
    cv = {}
    for categ in filter(lambda x: x!= 'genome', dvalues.keys()):
        cv[categ] = []
        for p in dvalues[categ]:
            allv.append([categ] + list(p))
            cv[categ].append(list(p))
    
    for p in sorted(allv, key=lambda x: (x[3], x[1], x[0]), reverse=True):
        fout.write( '\t'.join(p[1:3] + [p[0]] + [str(p[3])] + [str(p[6])]) + '\n')
    
    fout.close()
    logger.info('Saved combined pangenome informations (%s)'%fname)
    
def writeCombined(dvalues, orgs):
    '''
    Write down the table with the combined data
    '''
    fname = 'combined.tsv'
    fout = open(fname,'w')
    
    fout.write('# Metabolic reconstruction Vs. Phenomic activity\n')
    fout.write('# reactions = number of distinct reaction IDs\n')
    fout.write('# meanAV = mean AV\n')
    fout.write('# (-1, no activity available)\n')
    fout.write('\t'.join( ['', '', ''] + ['%s\t'%x for x in orgs] )
                    + '\n')
    fout.write('\t'.join( ['pathway', 'name', 'category'] +
                          ['reactions\tmeanAV' for x in orgs])
                    + '\n')
    
    # Order the pathway list for genomic content
    dallv = {}
    pnames = {}
    for categ in filter(lambda x: x!= 'genome', dvalues[orgs[0]].keys()):
        dallv[categ] = {}
        for org_id in orgs:
            for p in dvalues[org_id][categ]:
                dallv[categ][p[0]] = dallv[categ].get(p[0], [])
                dallv[categ][p[0]].append( (p[2],p[3]) )
                pnames[p[0]] = p[1]
    
    allv = []
    cv = {}
    for categ in filter(lambda x: x!= 'genome', dvalues[orgs[0]].keys()):
        cv[categ] = []
        for p in dallv[categ]:
            allv.append([categ, p, pnames[p]] + dallv[categ][p])
            cv[categ].append([p, pnames[p]] + dallv[categ][p])
    
    for p in sorted(allv, key=lambda x: [x[i][0] for i in range(3,len(orgs)+3)] + [x[1], x[0]], reverse=True):
        fout.write( '\t'.join(p[1:3] + [p[0]] +
                              ['%s\t%s'%(str(x[0]),str(x[1])) for x in p[3:]])
                    + '\n')
    
    fout.close()
    logger.info('Saved combined informations (%s)'%fname)

def dNet(project, allorgs=False, allpaths=False):
    '''
    Metabolic network reconstruction and analysis
    '''
    from ductape.common.utils import makeRoom
    from ductape.kegg.kegg import avoidedPaths
    
    kind = dSetKind(project)
    
    kegg = Kegg(project)
    
    # Check genomic and phenomic state
    proj = Project(project)
    proj.getProject()
    
    if proj.genome != 'map2kegg':
        logger.warning('Genome must be mapped to Kegg! (run dgenome start)')
        return True
    
    phenome = True
    if proj.phenome != 'map2kegg':
        logger.warning('Phenome is not mapped to Kegg, skipping that part')
        phenome = False
        
    biolog = Biolog(project)
    if biolog.atLeastOneNoParameter():
        logger.warning('Phenome parametrization has not yet been performed!')
        phenome = False
    
    logger.info('Saving overall metabolic network')
    aNet = getTotalNet(project)
    dapNet = {}
    for path in kegg.getMappedPathways():
        if path.path_id in avoidedPaths:continue
        dapNet[path.path_id] = getTotalNet(project, path.path_id)
        
    # Write
    npath = makeRoom('', 'metNet', 'KEGG')
    writeNet(aNet, npath, 'ALL.gml')
    for k,v in dapNet.iteritems():
        if ':' in k:
            k = k.split(':')[1]
        if allpaths:
            writeNet(v, npath, '%s.gml'%k)
        
    if proj.isPanGenome() and kind == 'pangenome' and not allorgs:
        orgs = ['conserved', 'variable']
        
        slen = 'metNet_pangenome_length.tsv'
        flen = open(slen,'w')
        flen.write('# Metabolic network length (number of distinct reaction IDs)\n')
        flen.write('\t'.join( ['network', 'name', 'overall'] + orgs) + '\n')
        
        sconn = 'metNet_pangenome_connected.tsv'
        fconn = open(sconn,'w')
        fconn.write('# Subnetworks (Connected components)\n')
        fconn.write('\t'.join( ['', ''] + ['Subnetworks'] +
                                              ['' for x in range(len(orgs))] +
                                              ['Subnetworks mean length'] +
                                              ['' for x in range(len(orgs))] +
                                              ['Subnetworks length std-dev'] +
                                              ['' for x in range(len(orgs))]) + '\n')
        fconn.write('\t'.join( ['network', 'name', 'overall'] + orgs +
                                                  ['overall'] + orgs + 
                                                  ['overall'] + orgs) + '\n')
        
        if phenome:
            sact = 'metNet_pangenome_activity.tsv'
            fact = open(sact,'w')
            fact.write('# Metabolic network activity\n')
            fact.write('\t'.join( ['network', 'name', 'category'] +
                                                 ['Avg Activity Difference'] +
                                                 ['Std Activity Difference']) + '\n')

        # Total network
        logger.info('Overall network stats')
        
        allr = [r for r in kegg.getMappedRPairsReact()]
        ecore = [r for r in kegg.getConservedRPairsReact()]
        edisp = [r for r in kegg.getVariableRPairsReact()]
        dpangenome = {'all': allr,'conserved':ecore, 'variable':edisp}
    
        oNet = {}
        for org in orgs:
            oNet[org] = getPanGenomeNet(project, dpangenome, org)
            npath = makeRoom('', 'metNet', org)
            writeNet(oNet[org], npath, '%s.gml'%org)
        
            
        flen.write('\t'.join( ['All', ''] +
                              [str(len(aNet.getDistinctReactions()))] +
                              [str(len(oNet[x].getDistinctReactions())) for x in orgs]) + '\n')
        
        fconn.write('\t'.join( ['All', ''] +
                              [str(aNet.getComponents())] +
                              [str(oNet[x].getComponents()) for x in orgs] +
                              [str(aNet.getComponentsMean())] +
                              [str(oNet[x].getComponentsMean()) for x in orgs] +
                              [str(aNet.getComponentsStd())] +
                              [str(oNet[x].getComponentsStd()) for x in orgs]) + '\n')
                              
        if phenome:
            for categ in biolog.getCategs(True):
                scateg = categ.category.replace(' ','_').replace('&','and')
                
                oNet = getPanGenomeNet(project,
                                       dpangenome,
                                       'all',
                                       category=categ.category)
                npath = makeRoom('', 'metNet', 'all', scateg)
                writeNet(oNet, npath, '%s_%s.gml'%('all', scateg))
                
                fact.write('\t'.join( ['All', '', scateg] +
                              [str(oNet.mean())] + [str(oNet.std())]) + '\n')
                              
        # Single paths
        logger.info('Single pathways stats')
        
        dPaths = {}
        
        for path in kegg.getMappedPathways():
            if path.path_id in avoidedPaths:continue
            logger.info('Pathway: %s // %s'%(path.path_id, path.name))
            
            if ':' in path.path_id:
                spath = path.path_id.split(':')[1]
            
            allr = [r for r in kegg.getMappedRPairsReact(path.path_id)]
            ecore = [r for r in kegg.getConservedRPairsReact(path.path_id)]
            edisp = [r for r in kegg.getVariableRPairsReact(path.path_id)]
            dpangenome = {'all': allr,'conserved':ecore, 'variable':edisp}
            
            paNet = getPanGenomeNet(project, dpangenome,
                                    'all', path_id=path.path_id)
            
            skip = False
            if len(paNet.getDistinctReactions()) == 0:
                logger.debug('Skipping reaction data on pathway: %s'%(path.path_id))
                skip = True
               
            oNet = {}
            for org_id in orgs:
                oNet[org_id] = getPanGenomeNet(project, dpangenome,
                                               org_id, path.path_id)
                if allpaths and not skip:
                    npath = makeRoom('', 'metNet', org_id)
                    writeNet(oNet[org_id], npath, '%s_%s.gml'%(org_id, spath))
            
            iAll = len(dapNet[path.path_id])
            iDisp = len(oNet['variable'].getDistinctReactions())
            iCore = len(oNet['conserved'].getDistinctReactions())
            iTotal = len(paNet.getDistinctReactions())
            
            if not skip:
                flen.write('\t'.join( [path.path_id, path.name] +
                                      [str(len(dapNet[path.path_id].getDistinctReactions()))] +
                                      [str(len(oNet[x].getDistinctReactions())) for x in orgs]) + '\n')
                
                fconn.write('\t'.join( [path.path_id, path.name] +
                                      [str(dapNet[path.path_id].getComponents())] +
                                      [str(oNet[x].getComponents()) for x in orgs] +
                                      [str(dapNet[path.path_id].getComponentsMean())] +
                                      [str(oNet[x].getComponentsMean()) for x in orgs] +
                                      [str(dapNet[path.path_id].getComponentsStd())] +
                                      [str(oNet[x].getComponentsStd()) for x in orgs]) + '\n')
            
            dAV = {}
            if phenome:
                path_co = [x.co_id for x in kegg.getPathComps(path.path_id)]
                for categ in biolog.getCategs(True):
                    wells = [w for w in biolog.getAllCoByCateg(categ.category)
                        if 'cpd:'+w.co_id in path_co]
                    scateg = categ.category.replace(' ','_').replace('&','and')
                    
                    if len(wells) == 0:
                        logger.debug('Skipping activity data on pathway: %s'
                                     %(path.path_id))
                        dAV[scateg] = -1
                        continue
                    
                    oNet = getPanGenomeNet(project,
                                     dpangenome,
                                     'all',
                                     path.path_id,
                                     categ.category)
                    
                    fAV = oNet.mean()
                    if math.isnan(fAV):
                        dAV[scateg] = -1
                    else:
                        dAV[scateg] = fAV
                    
                    if not oNet.hasNodesWeight():
                        logger.debug('Skipping activity data on pathway: %s (%s)'
                                    %(path.path_id, scateg))
                        continue
                    
                    if allpaths:
                        npath = makeRoom('', 'metNet', 'all', scateg)
                        writeNet(oNet, npath,
                                 '%s_%s_%s.gml'%('all', scateg, spath))
                    
                    fact.write('\t'.join( [path.path_id, path.name, scateg] +
                                  [str(oNet.mean())] + [str(oNet.std())]) + '\n')
                
            try:
                dtot = float(iDisp) / float(iTotal)
            except ZeroDivisionError:
                dtot = float('inf')
            try:
                ctot = float(iCore) / float(iTotal)
            except ZeroDivisionError:
                ctot = float('inf')
            try:
                dc = float(iDisp) / float(iCore)
            except ZeroDivisionError:
                dc = float('inf')
            
            dPaths['genome'] = dPaths.get('genome', [])
            dPaths['genome'].append((path.path_id, path.name, dtot, ctot, dc))
            for scateg in dAV:
                dPaths[scateg] = dPaths.get(scateg, [])
                dPaths[scateg].append((path.path_id, path.name, dtot, ctot,
                                       dc, dAV[scateg]))

    elif kind == 'single' or allorgs and not kind == 'mutants':
        logger.info('Single organisms networks')
        
        organism = Organism(project)
        
        orgs = [org.org_id for org in organism.getAll()]
        
        slen = 'metNet_length.tsv'
        flen = open(slen,'w')
        flen.write('# Metabolic network length (number of distinct reactions)\n')
        flen.write('\t'.join( ['network', 'name', 'overall'] + orgs) + '\n')
        
        sconn = 'metNet_connected.tsv'
        fconn = open(sconn,'w')
        fconn.write('# Subnetworks (Connected components)\n')
        fconn.write('\t'.join( ['', ''] + ['Subnetworks'] +
                                              ['' for x in range(len(orgs))] +
                                              ['Subnetworks mean length'] +
                                              ['' for x in range(len(orgs))] +
                                              ['Subnetworks length std-dev'] +
                                              ['' for x in range(len(orgs))]) + '\n')
        fconn.write('\t'.join( ['network', 'name', 'overall'] + orgs +
                                                  ['overall'] + orgs + 
                                                  ['overall'] + orgs) + '\n')
        
        if phenome:
            sact = 'metNet_activity.tsv'
            fact = open(sact,'w')
            fact.write('# Metabolic network activity\n')
            fact.write('\t'.join( ['', '', ''] + ['Mean Activity'] +
                                              ['' for x in range(len(orgs)-1)] +
                                              ['Activity std-dev'] +
                                              ['' for x in range(len(orgs)-1)])
                                            + '\n')
            fact.write('\t'.join( ['network', 'name', 'category'] +
                                                     orgs +
                                                     orgs) + '\n')
        
        # Total network
        logger.info('Overall network stats')
        
        oNet = {}
        for org_id in orgs:
            oNet[org_id] = getOrgNet(project, org_id)
            npath = makeRoom('', 'metNet', org_id)
            writeNet(oNet[org_id], npath, '%s.gml'%org_id)
            
        flen.write('\t'.join( ['All', ''] +
                              [str(len(aNet.getDistinctReactions()))] +
                              [str(len(oNet[x].getDistinctReactions())) for x in orgs])
                               + '\n')
        
        fconn.write('\t'.join( ['All', ''] +
                              [str(aNet.getComponents())] +
                              [str(oNet[x].getComponents()) for x in orgs] +
                              [str(aNet.getComponentsMean())] +
                              [str(oNet[x].getComponentsMean()) for x in orgs] +
                              [str(aNet.getComponentsStd())] +
                              [str(oNet[x].getComponentsStd()) for x in orgs])
                              + '\n')
        
        if phenome:
            for categ in biolog.getCategs(True):
                scateg = categ.category.replace(' ','_').replace('&','and')
                wells = [w for w in biolog.getAllCoByCateg(categ.category)]
                if len(wells) == 0:
                    logger.debug('Skipping activity data on pathway: %s'
                                 %(path.path_id))
                    continue
                
                oNet = {}
                for org_id in orgs:
                    oNet[org_id] = getOrgNet(project,
                                             org_id,
                                             category=categ.category)
                    npath = makeRoom('', 'metNet', org_id, scateg)
                    writeNet(oNet[org_id], npath, '%s_%s.gml'%(org_id, scateg))
                
                fact.write('\t'.join( ['All', '', scateg] +
                              [str(oNet[x].mean()) for x in orgs] +
                              [str(oNet[x].std()) for x in orgs]) + '\n')
        
        # Single paths
        logger.info('Single pathways stats')
        
        dPaths = {}
        for org_id in orgs:
            dPaths[org_id] = {}
        genome = {}
        
        for path in kegg.getMappedPathways():
            if path.path_id in avoidedPaths:continue
            logger.info('Pathway: %s // %s'%(path.path_id, path.name))
            
            if ':' in path.path_id:
                spath = path.path_id.split(':')[1]
                
            oNet = {}
            for org_id in orgs:
                oNet[org_id] = getOrgNet(project, org_id, path.path_id)
                
                dPaths[org_id]['genome'] = dPaths[org_id].get('genome', [])
                dPaths[org_id]['genome'].append([path.path_id, path.name,
                                                 len(oNet[org_id])])
                genome[org_id] = [path.path_id, path.name, len(oNet[org_id])]
                
                if len(oNet[org_id]) == 0:
                    logger.debug('Skipping reactions data on pathway: %s (%s)'
                                     %(path.path_id, org_id))
                    continue
                
                if allpaths:
                    npath = makeRoom('', 'metNet', org_id)
                    writeNet(oNet[org_id], npath, '%s_%s.gml'%(org_id, spath))
            
            skip = False
            if sum( [len(oNet[x]) for x in oNet] ) == 0:
                skip = True
                logger.debug('Skipping reactions data on pathway: %s'
                                     %(path.path_id))
            
            if not skip:
                flen.write('\t'.join( [path.path_id, path.name] +
                                  [str(len(dapNet[path.path_id].getDistinctReactions()))] +
                                  [str(len(oNet[x].getDistinctReactions())) for x in orgs])
                                  + '\n')
                
                fconn.write('\t'.join( [path.path_id, path.name] +
                                  [str(dapNet[path.path_id].getComponents())] +
                                  [str(oNet[x].getComponents()) for x in orgs] +
                                  [str(dapNet[path.path_id].getComponentsMean())] +
                                  [str(oNet[x].getComponentsMean()) for x in orgs] +
                                  [str(dapNet[path.path_id].getComponentsStd())] +
                                  [str(oNet[x].getComponentsStd()) for x in orgs]) + '\n')
            
            dAV = {}
            for org_id in orgs:
                dAV[org_id] = {}
                
            if phenome:
                path_co = [x.co_id for x in kegg.getPathComps(path.path_id)]
                for categ in biolog.getCategs(True):
                    wells = [w for w in biolog.getAllCoByCateg(categ.category)
                        if 'cpd:'+w.co_id in path_co]
                    
                    scateg = categ.category.replace(' ','_').replace('&','and')
                    if len(wells) == 0:
                        logger.debug('Skipping activity data on pathway: %s (%s)'
                                     %(path.path_id, scateg))
                        for org_id in orgs:
                            dAV[org_id][scateg] = -1
                        continue
                    
                    oNet = {}
                    for org_id in orgs:
                        oNet[org_id] = getOrgNet(project,
                                                 org_id,
                                                 path.path_id,
                                                 categ.category)
                        
                        fAV = oNet[org_id].mean()
                        if math.isnan(fAV):
                            dAV[org_id][scateg] = -1
                        else:
                            dAV[org_id][scateg] = fAV
                        
                        if allpaths:
                            if not oNet[org_id].hasNodesWeight():
                                logger.debug('Skipping activity data on pathway: %s (%s %s)'
                                     %(path.path_id, org_id, scateg))
                                continue
                            npath = makeRoom('', 'metNet', org_id, scateg)
                            writeNet(oNet[org_id], npath,
                                     '%s_%s_%s.gml'%(org_id, scateg, spath))
                    
                    check = set([oNet[x].hasNodesWeight() for x in oNet])
                    if len(check) == 1 and check.pop() == False:
                        logger.debug('Skipping activity data on pathway: %s (%s)'
                                     %(path.path_id, scateg))
                        continue
                    
                    fact.write('\t'.join( [path.path_id, path.name, scateg] +
                                  [str(oNet[x].mean()) for x in orgs] +
                                  [str(oNet[x].std()) for x in orgs]) + '\n')
                    
            for org_id in dPaths:
                for scateg in dAV[org_id]:
                    dPaths[org_id][scateg] = dPaths[org_id].get(scateg, [])
                    dPaths[org_id][scateg].append(genome[org_id] + [dAV[org_id][scateg]])
                    
    elif kind == 'mutants' or allorgs:
        logger.info('Mutants networks')
        
        organism = Organism(project)
        
        refs = [org.org_id
                    for org in organism.getAll()
                    if not organism.isMutant(org.org_id)]
        
        orgs = []
        ref_rpairs = {}
        
        for ref_id in refs:
            muts = [x for x in organism.getOrgMutants(ref_id)]
            ref_rpairs[ref_id] = kegg.getExclusiveRPairsReactMutants(ref_id, muts)
        
            orgs.append(ref_id)
            for m in muts:
                orgs.append(m)
        
        slen = 'metNet_length.tsv'
        flen = open(slen,'w')
        flen.write('# Metabolic network length (number of exclusive reactions w/ respect to the wild-type)\n')
        flen.write('\t'.join( ['network', 'name', 'overall'] + orgs) + '\n')
        
        sconn = 'metNet_connected.tsv'
        fconn = open(sconn,'w')
        fconn.write('# Subnetworks (Connected components)\n')
        fconn.write('\t'.join( ['', ''] + ['Subnetworks'] +
                                              ['' for x in range(len(orgs))] +
                                              ['Subnetworks mean length'] +
                                              ['' for x in range(len(orgs))] +
                                              ['Subnetworks length std-dev'] +
                                              ['' for x in range(len(orgs))]) + '\n')
        fconn.write('\t'.join( ['network', 'name', 'overall'] + orgs +
                                                  ['overall'] + orgs + 
                                                  ['overall'] + orgs) + '\n')
        
        if phenome:
            sact = 'metNet_activity.tsv'
            fact = open(sact,'w')
            fact.write('# Metabolic network activity\n')
            fact.write('\t'.join( ['', '', ''] + ['Mean Activity'] +
                                              ['' for x in range(len(orgs)-1)] +
                                              ['Activity std-dev'] +
                                              ['' for x in range(len(orgs)-1)])
                                            + '\n')
            fact.write('\t'.join( ['network', 'name', 'category'] +
                                                     orgs +
                                                     orgs) + '\n')
        
        # Total network
        logger.info('Overall network stats')
        
        oNet = {}
        
        for ref_id in refs:
            muts = [x for x in organism.getOrgMutants(ref_id)]
            
            oNet[ref_id] = getOrgNet(project, ref_id)
            npath = makeRoom('', 'metNet', ref_id)
            writeNet(oNet[ref_id], npath, '%s.gml'%ref_id)
            
            for mut_id in muts:
                oNet[mut_id] = getMutNet(project, mut_id,
                                         ref_rpairs[ref_id][mut_id].values())
                npath = makeRoom('', 'metNet', mut_id)
                writeNet(oNet[mut_id], npath, '%s.gml'%mut_id)
            
        flen.write('\t'.join( ['All', ''] +
                              [str(len(aNet.getDistinctReactions()))] +
                              [str(len(oNet[x].getDistinctReactions())) for x in orgs])
                              + '\n')
        
        fconn.write('\t'.join( ['All', ''] +
                              [str(aNet.getComponents())] +
                              [str(oNet[x].getComponents()) for x in orgs] +
                              [str(aNet.getComponentsMean())] +
                              [str(oNet[x].getComponentsMean()) for x in orgs] +
                              [str(aNet.getComponentsStd())] +
                              [str(oNet[x].getComponentsStd()) for x in orgs])
                              + '\n')
        
        if phenome:
            for categ in biolog.getCategs(True):
                scateg = categ.category.replace(' ','_').replace('&','and')
                wells = [w for w in biolog.getAllCoByCateg(categ.category)]
                if len(wells) == 0:
                    logger.debug('Skipping activity data on pathway: %s'
                                 %(path.path_id))
                    continue
                
                oNet = {}
                for ref_id in refs:
                    muts = [x for x in organism.getOrgMutants(ref_id)]
                    
                    oNet[ref_id] = getOrgNet(project, ref_id,
                                             category=categ.category)
                    npath = makeRoom('', 'metNet', ref_id, scateg)
                    writeNet(oNet[ref_id], npath, '%s_%s.gml'%(ref_id, scateg))
                    
                    for mut_id in muts:
                        oNet[mut_id] = getMutNet(project,
                                                 mut_id,
                                                 ref_rpairs[ref_id][mut_id].values(),
                                                 category=categ.category)
                        npath = makeRoom('', 'metNet', mut_id, scateg)
                        writeNet(oNet[mut_id], npath, '%s_%s.gml'%(mut_id, scateg))
                
                fact.write('\t'.join( ['All', '', scateg] +
                              [str(oNet[x].mean()) for x in orgs] +
                              [str(oNet[x].std()) for x in orgs]) + '\n')
        
        # Single paths
        logger.info('Single pathways stats')
        
        dPaths = {}
        for org_id in orgs:
            dPaths[org_id] = {}
        genome = {}
        
        for path in kegg.getMappedPathways():
            if path.path_id in avoidedPaths:continue
            logger.info('Pathway: %s // %s'%(path.path_id, path.name))
            
            if ':' in path.path_id:
                spath = path.path_id.split(':')[1]
                
            oNet = {}
        
            for ref_id in refs:
                muts = [x for x in organism.getOrgMutants(ref_id)]
                ref_rpairs = kegg.getExclusiveRPairsReactMutants(ref_id, muts,
                                                                 path.path_id)
            
                oNet[ref_id] = getOrgNet(project, ref_id, path.path_id)
                
                dPaths[ref_id]['genome'] = dPaths[ref_id].get('genome', [])
                dPaths[ref_id]['genome'].append([path.path_id, path.name,
                                                 len(oNet[ref_id])])
                genome[ref_id] = [path.path_id, path.name, len(oNet[ref_id])]
                
                skip = False
                if len(oNet[ref_id]) == 0:
                    logger.debug('Skipping reactions data on pathway: %s (%s)'
                                     %(path.path_id, ref_id))
                    skip = True
                
                if allpaths and not skip:
                    npath = makeRoom('', 'metNet', ref_id)
                    writeNet(oNet[ref_id], npath, '%s_%s.gml'%(ref_id, spath))
            
                for mut_id in muts:
                    oNet[mut_id] = getMutNet(project,
                                             mut_id,
                                             ref_rpairs[mut_id].values(),
                                             path.path_id)
                
                    dPaths[mut_id]['genome'] = dPaths[mut_id].get('genome', [])
                    dPaths[mut_id]['genome'].append([path.path_id, path.name,
                                                     len(oNet[mut_id])])
                    genome[mut_id] = [path.path_id, path.name, len(oNet[mut_id])]
                
                    if len(oNet[mut_id]) == 0:
                        logger.debug('Skipping reactions data on pathway: %s (%s)'
                                         %(path.path_id, mut_id))
                        continue
                    
                    if allpaths and not skip:
                        npath = makeRoom('', 'metNet', mut_id)
                        writeNet(oNet[mut_id], npath, '%s_%s.gml'%(mut_id, spath))
            
            skip = False
            if sum( [len(oNet[x]) for x in oNet] ) == 0:
                skip = True
                logger.debug('Skipping reactions data on pathway: %s'
                                     %(path.path_id))
            
            if not skip:
                flen.write('\t'.join( [path.path_id, path.name] +
                                  [str(len(dapNet[path.path_id].getDistinctReactions()))] +
                                  [str(len(oNet[x].getDistinctReactions())) for x in orgs])
                                  + '\n')
                
                fconn.write('\t'.join( [path.path_id, path.name] +
                                  [str(dapNet[path.path_id].getComponents())] +
                                  [str(oNet[x].getComponents()) for x in orgs] +
                                  [str(dapNet[path.path_id].getComponentsMean())] +
                                  [str(oNet[x].getComponentsMean()) for x in orgs] +
                                  [str(dapNet[path.path_id].getComponentsStd())] +
                                  [str(oNet[x].getComponentsStd()) for x in orgs])
                                  + '\n')
            
            dAV = {}
            for org_id in orgs:
                dAV[org_id] = {}
            
            if phenome:
                path_co = [x.co_id for x in kegg.getPathComps(path.path_id)]
                for categ in biolog.getCategs(True):
                    wells = [w for w in biolog.getAllCoByCateg(categ.category)
                        if 'cpd:'+w.co_id in path_co]
                    
                    scateg = categ.category.replace(' ','_').replace('&','and')
                    if len(wells) == 0:
                        logger.debug('Skipping activity data on pathway: %s (%s)'
                                     %(path.path_id, scateg))
                        for org_id in orgs:
                            dAV[org_id][scateg] = -1
                        continue
                    
                    oNet = {}
                    
                    for ref_id in refs:
                        muts = [x for x in organism.getOrgMutants(ref_id)]
                        ref_rpairs = kegg.getExclusiveRPairsReactMutants(ref_id, muts,
                                                                         path.path_id)
                    
                        oNet[ref_id] = getOrgNet(project,
                                                 ref_id,
                                                 path.path_id,
                                                 categ.category)
                        
                        fAV = oNet[ref_id].mean()
                        if math.isnan(fAV):
                            dAV[ref_id][scateg] = -1
                        else:
                            dAV[ref_id][scateg] = fAV
                        
                        if allpaths:
                            skip = False
                            if not oNet[ref_id].hasNodesWeight():
                                logger.debug('Skipping activity data on pathway: %s (%s %s)'
                                     %(path.path_id, ref_id, scateg))
                                skip = True
                            
                            if not skip:
                                npath = makeRoom('', 'metNet', ref_id, scateg)
                                writeNet(oNet[ref_id], npath,
                                         '%s_%s_%s.gml'%(ref_id, scateg, spath))
                        
                        for mut_id in muts:
                            oNet[mut_id] = getMutNet(project,
                                                 mut_id,
                                                 ref_rpairs[mut_id].values(),
                                                 path.path_id,
                                                 categ.category)
                            
                            fAV = oNet[mut_id].mean()
                            if math.isnan(fAV):
                                dAV[mut_id][scateg] = -1
                            else:
                                dAV[mut_id][scateg] = fAV
                        
                            if allpaths:
                                if not oNet[mut_id].hasNodesWeight():
                                    logger.debug('Skipping activity data on pathway: %s (%s %s)'
                                         %(path.path_id, mut_id, scateg))
                                    continue
                                npath = makeRoom('', 'metNet', mut_id, scateg)
                                writeNet(oNet[mut_id], npath,
                                         '%s_%s_%s.gml'%(mut_id, scateg, spath))
                    
                    check = set([oNet[x].hasNodesWeight() for x in oNet])
                    if len(check) == 1 and check.pop() == False:
                        logger.debug('Skipping activity data on pathway: %s (%s)'
                                     %(path.path_id, scateg))
                        continue
                    
                    if not skip:
                        fact.write('\t'.join( [path.path_id, path.name, scateg] +
                                      [str(oNet[x].mean()) for x in orgs] +
                                      [str(oNet[x].std()) for x in orgs]) + '\n')
            
            for org_id in dPaths:
                for scateg in dAV[org_id]:
                    dPaths[org_id][scateg] = dPaths[org_id].get(scateg, [])
                    dPaths[org_id][scateg].append(genome[org_id] + [dAV[org_id][scateg]])
    
    try:
        flen.close()
        fconn.close()
        if phenome:
            fact.close()
    except:pass
    
    logger.info('Metabolic network length stats saved to %s'%slen)
    logger.info('Metabolic subnetworks stats saved to %s'%sconn)
    if phenome:
        logger.info('Metabolic network activity stats saved to %s'%sact)
    
    logger.info('Saving combined metrics and stats')
    if proj.isPanGenome() and kind == 'pangenome' and not allorgs:
        writeCombinedPanGenome(dPaths)
    elif kind == 'single' or allorgs or kind == 'mutants':
        writeCombined(dPaths, orgs)
    return True

def getCombinedMatrix(phenome, genome, matrix, pthresh, gthresh):
    '''
    Given a phenomic data and genomic data vectors and a data dictionary
    return a plottable matrix w/ compounds on y-axis and w/ pathways on x-axis
    
    Labels are returned as well
    
    thresholds are inclusive
    '''
    import numpy as np

    phenome = sorted(filter(lambda x: x[3] >= pthresh, phenome),
                     key=lambda x: x[3], reverse=True)
    pnames = [(x[1]+' '+x[0], x[2]) for x in phenome]
    
    genome = sorted(filter(lambda x: x[2] >= 0 , genome),
                    key=lambda x: x[2])
    
    # Remove those pathways with no compound mapped
    toRemove = []
    for j in genome:
        p = j[0]
        
        bOne = False
        for i in phenome:
            cid = i[0]
            scateg = i[1]

            if p in matrix[cid][scateg]:
                bOne = True
                break
        
        if not bOne:
            toRemove.append(j)
    for i in toRemove:
        genome.remove(i)
    
    gnames = [(x[0], x[1]) for x in genome]

    matr=[]
    for i in phenome:
        vec = []
        
        cid = i[0]
        scateg = i[1]
        
        for j in genome:
            p = j[0]
            if p in matrix[cid][scateg]:
                vec.append(matrix[cid][scateg][p])
            else:
                vec.append(np.nan)
        
        matr.append(vec)
    
    return matr, pnames, gnames

def writeCombinedMatrix(fhandle, matr, pnames, gnames):
    '''
    Given a file handle, a matrix and relative labels print out a combined matrix
    '''
    fhandle.write('\t'.join([' '] + [x[0] + ' ' + x[1] for x in gnames]))
    fhandle.write('\n')
    
    for i in range(len(matr)):
        v = matr[i]
        fhandle.write('\t'.join([pnames[i][0] + ' ' + pnames[i][1]]
                                + [str(x) for x in v]))
        fhandle.write('\n')
        
def plotCombinedMatrix(fname, matr, pnames, gnames, cmap=None,
                       vmin=0, vmax=None,
                       ylabel='Phenotypic activity',
                       xlabel='Pathways',
                       title='Combined genome/phenome'):
    '''
    Given a file name, a matrix and relative labels plot out a combined heatmap
    
    also the colormap and min and max values have to be provided
    '''
    if len(matr) == 0:
        logger.warning('No data available for a combined plot')
        return
    
    if len(gnames) < 50:
        w = 15
    else:
        w = len(gnames)/3
    
    if len(pnames) < 50:
        h = 15
    else:
        h = len(pnames)/3
        
    fig = plt.figure(figsize=(w,h))
    ax = fig.add_subplot(111)
    
    if not cmap:
        cmap=cm.RdYlGn
    cmap.set_bad('gray',.4)
    
    cb = ax.imshow(matr, aspect='auto', interpolation='none',
              cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_yticks(np.arange(len(pnames)))
    ax.set_yticklabels([x[0] for x in pnames], size=9)
    ax.set_xticks(np.arange(len(gnames)))
    ax.set_xticklabels([x[0] for x in gnames], size=9,
                  rotation=90)
    
    if w > 15 or h > 15:
        plt.colorbar(cb, shrink=0.5, aspect=30)
    else:
        plt.colorbar(cb)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    ax.set_title(title)
    
    plt.tight_layout()
    
    plt.savefig(fname)
    
def getCombinations(matrix, phenome, genome, pthresh, gthresh):
    '''
    Generator to combination of compounds and pathways data
    '''
    phenome = sorted(filter(lambda x: x[3] >= pthresh, phenome),
                     key=lambda x: x[3], reverse=True)
    
    genome = sorted(filter(lambda x: x[2] >= 0 , genome),
                    key=lambda x: x[2])
    
    # Remove those pathways with no compound mapped
    toRemove = []
    for j in genome:
        p = j[0]
        
        bOne = False
        for i in phenome:
            cid = i[0]
            scateg = i[1]

            if p in matrix[cid][scateg]:
                bOne = True
                break
        
        if not bOne:
            toRemove.append(j)
    for i in toRemove:
        genome.remove(i)
    
    out = []
    for i in phenome:
        cid = i[0]
        scateg = i[1]
        cname = i[2]
        pval = i[3]
        
        for j in genome:
            p = j[0]
            pname = j[1]
            gval = j[2]
            
            if p in matrix[cid][scateg]:
                out.append((scateg, cid, cname, p, pname, pval, gval))
    
    for o in sorted(out, key=lambda x: (x[5], x[6]), reverse=True):
        yield o

def dCombine(project, allorgs=False, pthresh=5, doPrint=True):
    '''
    Prepare a table/heatmap focused on compound activity/genetic content
    '''
    from ductape.kegg.kegg import avoidedPaths
    from itertools import combinations
    
    kind = dSetKind(project)
    
    kegg = Kegg(project)
    
    # Check genomic and phenomic state
    proj = Project(project)
    proj.getProject()
    
    if proj.genome != 'map2kegg':
        logger.warning('Genome must be mapped to Kegg! (run dgenome start)')
        return True
    
    if proj.phenome != 'map2kegg':
        logger.warning('Phenome must be mapped to Kegg!')
        return True
        
    biolog = Biolog(project)
    if biolog.atLeastOneNoParameter():
        logger.warning('Phenome parametrization has not yet been performed!')
        return True
    
    if proj.isPanGenome() and kind == 'pangenome' and not allorgs:
        logger.info('Analyzing combined data from a pangenome')
           
        matr = {}
        cos = []
        gens = []
        
        # Start from the pathways with at least one reaction mapped
        paths = set([x.path_id for x in kegg.getMappedPathways()])
        for p in avoidedPaths:
            try:
                paths.remove(p)
            except:pass
           
        # Cycle through category/compound to build the matrix
        for categ in biolog.getCategs(True):
            category = categ.category
            scateg = categ.category.replace(' ','_').replace('&','and')

            corg = {}
            
            wells = [w for w in biolog.getAllCoByCateg(category)]
            
            for well in wells:
                co_id = 'cpd:' + well.co_id
                
                acts = [x.avgact 
                        for x in
                        biolog.getAvgActivityEachOrg(well.plate_id,
                                                     well.well_id)]
                if len(acts) <= 1:
                    continue
                
                diffs = []
                for a, a1 in combinations(acts, 2):
                    diffs.append( abs(a - a1) )
                    
                avgdiff = np.array(diffs).mean()
                
                # Some co_ids are present more than once
                if co_id not in corg:
                    corg[co_id] = []
                corg[co_id].append(avgdiff)
            
            toremove = set()
            for k, v in corg.iteritems():
                mean = np.array(v).mean()
                if not math.isnan(float(mean)):
                    corg[k] = mean
                else:
                    toremove.add(k)
            for k in toremove:
                del corg[k]
                
            # Insert the compounds inside the matrix
            for co_id in corg:
                # Pathways mapped to this co_id
                # only those with reactions will be used
                for p in kegg.getCompPaths(co_id):
                    if p.path_id not in paths:continue
                    
                    matr[co_id] = matr.get(co_id, {})
                    matr[co_id][scateg] = matr[co_id].get(scateg, {})
                    matr[co_id][scateg][p.path_id] = matr[co_id][scateg].get(
                                                         p.path_id, corg[co_id])
                
                # Add the phenotypic variability
                if co_id in matr and scateg in matr[co_id]:
                    cname = kegg.getCompound(co_id).name
                    cos.append((co_id, scateg, cname, corg[co_id]))
        
        # Get the genetic variability
        for p in paths:
            allr = [r for r in kegg.getMappedRPairsReact(p)]
            ecore = [r for r in kegg.getConservedRPairsReact(p)]
            edisp = [r for r in kegg.getVariableRPairsReact(p)]
            dpangenome = {'all': allr,'conserved':ecore, 'variable':edisp}
                
            totNet = len(getPanGenomeNet(project,
                                     dpangenome, 'all',
                                     path_id=p).getDistinctReactions())
            dispNet = len(getPanGenomeNet(project,
                                      dpangenome, 'variable',
                                      path_id=p).getDistinctReactions())
            
            pname = kegg.getPathway(p).name
            try:
                gens.append((p, pname, float(dispNet)/float(totNet)))
            except:pass
        
        matr_all, phen_all, gen_all = getCombinedMatrix(cos, gens, matr,
                                                        0, 0)
            
        # Write the whole matrix
        fname = 'combined_matrix_full.tsv'
        fout = open(fname, 'w')
        fout.write('# Combined matrix for the pangenome\n')
        fout.write('# Each cell contains the average difference on the AV '+
                   ' between each organism of the pangenome\n')
        fout.write('#  Compounds are sorted by mean diffAV, '+
                   ' pathways by genomic variability '+
                   '(number of variable reaction IDs / '+
                   'numer of distinct reaction IDs)\n')
        writeCombinedMatrix(fout, matr_all, phen_all, gen_all)
        
        logger.info('Saved overall combined pangenome matrix (%s)'%fname)
        
        # Reduced matrix
        matr_red, phen, gen = getCombinedMatrix(cos, gens, matr,
                                            pthresh, 0.0000001)
        
        # Write the matrix
        fname = 'combined_matrix.tsv'
        fout = open(fname, 'w')
        fout.write('# Combined matrix for the pangenome\n')
        fout.write('# Each cell contains the average difference on the AV '+
                   ' between each organism of the pangenome\n')
        fout.write('#  Compounds are sorted by mean diffAV, '+
                   ' pathways by genomic variability '+
                   '(number of variable reaction IDs / '+
                   'numer of distinct reaction IDs)\n')
        fout.write('#  mean diffAV threshold: %f\n'%pthresh)
        fout.write('#  genomic variability threshold: > 0\n')
        writeCombinedMatrix(fout, matr_red, phen, gen)
        
        logger.info('Saved reduced combined pangenome matrix (%s)'%fname)
        
        # Plot!
        fname = 'combined_variability.png'
        plotCombinedMatrix(fname, matr_red, phen, gen, cmap=cm.Purples,
                       vmax=biolog.getMaxActivity(),
                       ylabel='Phenotypic activity (mean diffAV)',
                       title='Combined genome/phenome variability')
        
        logger.info('Saved combined pangenome plot (%s)'%fname)
        
        # Print relevant combinations
        logger.info('Relevant combined data')
        
        header = '\t'.join( ['category', 'co_id', 'name', 
                             'path_id', 'name',
                             'mean diffAV', 'genomic variability'] )
        if doPrint:
            print(header)
        else:
            logger.info(header)
        
        for scateg, cid, cname, p, pname, pval, gval in getCombinations(matr,
                                                                    cos,
                                                                    gens,
                                                                    pthresh,
                                                                    0.0000001):
            line = '\t'.join( [str(x)
                               for x in [scateg, cid, cname,
                                         p, pname, pval, gval]] )
            
            if doPrint:
                print(line)
            else:
                logger.info(line)
        
    elif kind == 'single' or allorgs and not kind == 'mutants':
        logger.info('Analyzing combined data for each organisms')
        
        organism = Organism(project)
        
        orgs = [org.org_id for org in organism.getAll()]
        
        for org_id in orgs:
            # Start from the pathways with at least one reaction mapped
            paths = set([x.path_id for x in kegg.getMappedPathways(org_id)])
            for p in avoidedPaths:
                try:
                    paths.remove(p)
                except:pass
            
            matr = {}
            cos = []
            gens = []
            # Cycle through category/compound to build the matrix
            for categ in biolog.getCategs(True):
                category = categ.category
                scateg = categ.category.replace(' ','_').replace('&','and')
    
                corg = {}
        
                # Filter by path?
                wells = [w for w in biolog.getAllCoByCateg(category)]
                for well in wells:
                    co_id = 'cpd:' + well.co_id
                    
                    act = biolog.getAvgActivity(well.plate_id, well.well_id,
                                                org_id)
                    if act is not None:
                        # Some co_ids are present more than once
                        if co_id not in corg:
                            corg[co_id] = []
                        corg[co_id].append(act)
            
                toremove = set()
                for k, v in corg.iteritems():
                    mean = np.array(v).mean()
                    if not math.isnan(float(mean)):
                        corg[k] = mean
                    else:
                        toremove.add(k)
                for k in toremove:
                    del corg[k]
                
                # Insert the compounds inside the matrix
                for co_id in corg:
                    # Pathways mapped to this co_id
                    # only those with reactions will be used
                    for p in kegg.getCompPaths(co_id):
                        if p.path_id not in paths:continue
                        
                        matr[co_id] = matr.get(co_id, {})
                        matr[co_id][scateg] = matr[co_id].get(scateg, {})
                        matr[co_id][scateg][p.path_id] = matr[co_id][scateg].get(
                                                             p.path_id, corg[co_id])
                    
                    # Add the phenotypic variability
                    if co_id in matr and scateg in matr[co_id]:
                        cname = kegg.getCompound(co_id).name
                        cos.append((co_id, scateg, cname, corg[co_id]))
            
            # Get the genetic content
            for p in paths:
                pname = kegg.getPathway(p).name
                gens.append((p, pname,
                     len(getOrgNet(project, org_id, path_id=p).getDistinctReactions())))
                
            matr_all, phen_all, gen_all = getCombinedMatrix(cos, gens, matr,
                                                        0, 0)
            
            # Write the whole matrix
            fname = 'combined_matrix_%s.tsv'%org_id
            fout = open(fname, 'w')
            fout.write('# Combined matrix for %s\n'%org_id)
            fout.write('# Each cell contains the AV '+
                       ' for each compound\n')
            fout.write('#  Compounds are sorted by AV, '+
                       ' pathways by genomic content '+
                       '(numer of distinct reaction IDs)\n')
            writeCombinedMatrix(fout, matr_all, phen_all, gen_all)
            
            logger.info('Saved overall combined matrix (%s)'%fname)
            
            # Plot!
            fname = 'combined_%s.png'%org_id
            plotCombinedMatrix(fname, matr_all, phen_all, gen_all,
                           vmax=biolog.getMaxActivity())
            
            logger.info('Saved combined genome/phenome plot (%s)'%fname)
        
    elif kind == 'mutants' or allorgs:
        logger.info('Analyzing combined data for each mutant w/r/t the wild-type')
    
        organism = Organism(project)
        
        refs = [org.org_id
                    for org in organism.getAll()
                    if not organism.isMutant(org.org_id)]
        
        for ref_id in refs:
            muts = [x for x in organism.getOrgMutants(ref_id)]
            
            for mut_id in muts:
                # Start from the pathways with at least one reaction mapped
                paths = set([x.path_id for x in kegg.getMappedPathways(mut_id)])
                for p in avoidedPaths:
                    try:
                        paths.remove(p)
                    except:pass
                
                matr = {}
                cos = []
                gens = []
                # Cycle through category/compound to build the matrix
                for categ in biolog.getCategs(True):
                    category = categ.category
                    scateg = categ.category.replace(' ','_').replace('&','and')

                    corg = {}
        
                    # Filter by path?
                    wells = [w for w in biolog.getAllCoByCateg(category)]
                    for well in wells:
                        co_id = 'cpd:' + well.co_id
                        refact = biolog.getAvgActivity(well.plate_id, well.well_id,
                                                       ref_id)
                        act = biolog.getAvgActivity(well.plate_id, well.well_id,
                                                    mut_id)
                        if act is not None and refact is not None:
                            # Some co_ids are present more than once
                            if co_id not in corg:
                                corg[co_id] = []
                            corg[co_id].append(refact-act)
                
                    toremove = set()
                    for k, v in corg.iteritems():
                        mean = np.array(v).mean()
                        if not math.isnan(float(mean)):
                            corg[k] = mean
                        else:
                            toremove.add(k)
                    for k in toremove:
                        del corg[k]
                    
                    # Insert the compounds inside the matrix
                    for co_id in corg:
                        # Pathways mapped to this co_id
                        # only those with reactions will be used
                        for p in kegg.getCompPaths(co_id):
                            if p.path_id not in paths:continue
                            
                            matr[co_id] = matr.get(co_id, {})
                            matr[co_id][scateg] = matr[co_id].get(scateg, {})
                            matr[co_id][scateg][p.path_id] = matr[co_id][scateg].get(
                                                                 p.path_id, corg[co_id])
                        
                        # Add the phenotypic variability
                        if co_id in matr and scateg in matr[co_id]:
                            cname = kegg.getCompound(co_id).name
                            cos.append((co_id, scateg, cname, corg[co_id]))
                
                # Get the genetic content
                for p in paths:
                    pname = kegg.getPathway(p).name
                    gens.append((p, pname,
                         len(getOrgNet(project, mut_id,
                                       path_id=p).getDistinctReactions())))
                    
                matr_all, phen_all, gen_all = getCombinedMatrix(cos, gens, matr,
                                                        0, 0)
            
                # Write the whole matrix
                fname = 'combined_matrix_full_%s.tsv'%mut_id
                fout = open(fname, 'w')
                fout.write('# Combined matrix for the %s mutant (WT %s)\n'%(mut_id,
                                                                            ref_id))
                fout.write('# Each cell contains the difference on the AV '+
                           ' between the mutant and the wild-type\n')
                fout.write('#  Compounds are sorted by diffAV, '+
                           ' pathways by genomic content '+
                           '(numer of distinct mutated reaction IDs)\n')
                writeCombinedMatrix(fout, matr_all, phen_all, gen_all)
                
                logger.info('Saved overall combined genome/phenome matrix (%s)'%fname)
                
                # Reduced matrix
                matr_comb, phen, gen = getCombinedMatrix(cos, gens, matr,
                                                    pthresh, 0.0000001)
                
                # Write the matrix
                fname = 'combined_matrix_%s.tsv'%mut_id
                fout = open(fname, 'w')
                fout.write('# Combined matrix for the %s mutant (WT %s)\n'%(mut_id,
                                                                            ref_id))
                fout.write('# Each cell contains the difference on the AV '+
                           ' between the mutant and the wild-type\n')
                fout.write('#  Compounds are sorted by diffAV, '+
                           ' pathways by genomic content '+
                           '(numer of distinct mutated reaction IDs)\n')
                fout.write('#  diffAV threshold: %f\n'%pthresh)
                writeCombinedMatrix(fout, matr_comb, phen, gen)
                
                logger.info('Saved reduced combined pangenome matrix (%s)'%fname)
                
                # Plot!
                fname = 'combined_%s.png'%mut_id
                plotCombinedMatrix(fname, matr_comb, phen, gen, cmap=cm.PuOr,
                               vmin=-biolog.getMaxActivity(),
                               vmax=biolog.getMaxActivity(),
                               xlabel='Pathways containing the mutated reactions',
                               ylabel='Phenotypic variability w/r/t wild-type (diffAV)')
                
                logger.info('Saved combined genome/phenome plot (%s)'%fname)
                
                # Print relevant combinations
                logger.info('Relevant combined data')
                
                header = '\t'.join( ['category', 'co_id', 'name', 
                                        'path_id', 'name',
                                        'diffAV',
                                        'distinct mutated reaction IDs'] )
                if doPrint:
                    print(header)
                else:
                    logger.info(header)
                
                for scateg, cid, cname, p, pname, pval, gval in getCombinations(matr,
                                                                        cos,
                                                                        gens,
                                                                        pthresh,
                                                                        0.0000001):
                    line = '\t'.join( [str(x)
                                       for x in [scateg, cid, cname,
                                                 p, pname, pval, gval]] )
                    
                    if doPrint:
                        print(line)
                    else:
                        logger.info(line)
    
    return True

def getPlatesOrder(project):
    from ductape.storage.SQLite.database import Biolog
    
    return [x.plate_id for x in Biolog(project).getPlates()]

def getOrder(project, plates=None):
    '''
    Generator of plate/well IDs
    If plates is provided as a list of plates IDs, only those plates are used
    '''
    from ductape.storage.SQLite.database import Biolog
    
    if not plates:
        plates = getPlatesOrder(project)
        
    for pid in plates:
        for wid in Biolog(project).getPlateWells(pid):
            yield (pid, wid)
 
def dSetKind(project):
    '''
    Set the kind of genomic project and return its value
    '''
    proj = Project(project)
    proj.getProject()
    org = Organism(project)
    if org.howManyMutants() > 0:
        logger.info('%d mutants are present'%org.howManyMutants())
        if proj.kind != 'mutants':
            proj.setKind('mutants')
        return 'mutants'
    
    elif org.howMany() == 1:
        logger.info('Just one organism is present')
        if proj.kind != 'single':
            proj.setKind('single')
        return 'single'
    
    elif org.howMany() == 0:
        logger.info('No organisms are present yet')
        return None
    
    else:
        logger.info('%d organisms are present'%org.howMany())
        if proj.kind != 'pangenome':
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
        
        return steps

def mergeKegg(project):
    '''
    Returns some info on new kegg annotations extracted from orthologs
    
    set with (prot_id, ko_id)
    set with group_id having re-annotated proteins
    set with group_id having multiple ko ids
    '''
    genome = Genome(project)
    pankegg = genome.getPanGenomeKOs()
    
    merged = set()
    mergedg = set()
    
    multiple = set()
    
    for group_id, pko in pankegg.iteritems():
        allko = set()
        missing = 0
        for prot_id, kos in pko.iteritems():
            if kos is None:
                missing += 1
                continue
            for ko in kos:
                allko.add(ko)
        
        # How many distinct IDs do we have
        if len(allko) == 1 and missing >= 1:
            mergedg.add(group_id)
            # Single-annotations missing
            for prot_id, kos in pko.iteritems():
                if kos is None:
                    for ko in allko:
                        merged.add((prot_id, ko))
                        
        elif len(allko) > 1:
            # Multiple annotations!
            if missing == 0:
                # Only duplicates
                multiple.add(group_id)
                for prot_id, kos in pko.iteritems():
                    if len(kos) != len(allko):
                        mergedg.add(group_id)
                        
                        for ko in allko:
                            if ko not in pko[prot_id]:
                                merged.add((prot_id, ko))
                                     
            else:
                # Missing (and?) multiple
                mergedg.add(group_id)
                multiple.add(group_id)
                
                for prot_id, kos in pko.iteritems():
                    if kos is None:
                        for ko in allko:
                            merged.add((prot_id, ko))
                    else:
                        if len(kos) != len(allko):
                            for ko in allko:
                                if ko not in pko[prot_id]:
                                    merged.add((prot_id, ko))
                                    
    return merged, mergedg, multiple

def getPathsReacts(project):
    from ductape.kegg.kegg import avoidedPaths
    kegg = Kegg(project)
    # Get the pathway - reaction links
    paths = {}
    for pR in kegg.getPathReacts():
        if pR.path_id in avoidedPaths:
            continue
        if pR.path_id not in paths:
            paths[pR.path_id] = []
        paths[pR.path_id].append(pR.re_id)
        
    return paths

def getPathsComps(project):
    from ductape.kegg.kegg import avoidedPaths
    kegg = Kegg(project)
    # Get the pathway - compounds links
    paths = {}
    for pR in kegg.getPathComps():
        if pR.path_id in avoidedPaths:
            continue
        if pR.path_id not in paths:
            paths[pR.path_id] = []
        paths[pR.path_id].append(pR.co_id)
        
    return paths

def getExclusiveReactions(project, orgs=set()):
    '''
    Given a bunch of organisms, get two dicts
    mreacts: org_id -->  set(re_id, ...) (common)
    ereacts: org_id -->  set(re_id, ...) (exclusive)
    '''
    from ductape.storage.SQLite.database import Kegg
    
    kegg = Kegg(project)
    
    # Get the exclusive reactions
    ereacts = kegg.getExclusiveReactions(orgs)
    # Get the reactions of each organisms
    mreacts = {}
    for org_id in orgs:
        mreacts[org_id] = set()
        for oR in kegg.getOrgReact(org_id):
            mreacts[org_id].add(oR.re_id)
    
    # Remove the exclusives
    for org_id in orgs:
        mreacts[org_id].difference_update(ereacts[org_id])
        
    return mreacts, ereacts

def getExclusiveReactionsMutants(project, ref_id, muts=set()):
    '''
    Given a bunch of organisms, get two dicts
    mreacts: org_id -->  set(re_id, ...) (wild-type)
    mix: org_id --> set(re_id, ...) (mutated but also in wild-type)
    ereacts: org_id -->  set(re_id, ...) (exclusive mutated)
    '''
    from ductape.storage.SQLite.database import Kegg
    
    kegg = Kegg(project)
    
    # Get the exclusive reactions
    ereacts = kegg.getExclusiveReactionsMutants(ref_id, muts)    
    # Get the reactions of each organisms
    mreacts = {}
    mix = {}
    for mut_id in muts:
        # Reference genome reactions
        mut = set()
        for oR in kegg.getOrgReact(mut_id):
            mut.add(oR.re_id)
        
        mreacts[mut_id] = set()
        mix[mut_id] = set()
        for oR in kegg.getReferenceReact(mut_id, ref_id):
            if oR.re_id in mut:
                mix[mut_id].add(oR.re_id)
            else:
                mreacts[mut_id].add(oR.re_id)
        
    return mreacts, mix, ereacts

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

def createLegend(kind, project):
    '''
    Create a color scheme legend
    '''
    # TODO: a more centralized color scheme
    fig = plt.figure()
    fname = 'legend.png'
    
    # Get the maximum activity now present
    maxAct = Biolog(project).getMaxActivity()
    if maxAct is None:
        maxAct = 9
    
    rmatrix = np.outer(np.arange(0.33,1,0.01),np.ones(7))
    if kind == 'mutants' or 'singlediff':
        cmatrix = np.outer(np.arange(-maxAct,maxAct,0.1),np.ones(7))
    elif kind == 'pangenome':
        cmatrix = np.outer(np.arange(0.33,1,0.01),np.ones(7))
    else:
        cmatrix = np.outer(np.arange(0,maxAct,0.1),np.ones(7))
         
    if kind == 'pangenome':
        ax = fig.add_subplot(131, axisbg='b')
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Conserved')
        
        ax = fig.add_subplot(132)
        ax.imshow(rmatrix, cmap=cm.autumn, vmin=0, vmax=1)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Variable')
        
        ax = fig.add_subplot(133)
        ax.imshow(cmatrix, cmap=cm.Purples, vmin=0, vmax=maxAct)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Mean activity diff.')
        
        fig.savefig(fname)
            
    elif kind == 'single':
        ax = fig.add_subplot(121, axisbg='b')
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Reactions')
        
        ax = fig.add_subplot(122)
        ax.imshow(cmatrix, cmap=cm.RdYlGn, vmin=0, vmax=maxAct)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Phenomic activity')
        
        fig.savefig(fname)
        
    elif kind == 'singlediff':
        ax = fig.add_subplot(121, axisbg='b')
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Reactions')
        
        ax = fig.add_subplot(122)
        ax.imshow(cmatrix, cmap=cm.PuOr, vmin=-maxAct, vmax=maxAct)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Delta-activity')
        
        fig.savefig(fname)

    elif kind == 'mutants':
        ax = fig.add_subplot(141, axisbg='b')
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Wild-type')
        
        ax = fig.add_subplot(142, axisbg='g')
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Wild-type and Mutated')
        
        ax = fig.add_subplot(143, axisbg=pltcls.cnames['yellow'])
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])
        ax.axes.get_xaxis().set_ticks([])
        ax.set_title('Mutated')
        
        ax = fig.add_subplot(144)
        ax.imshow(cmatrix, cmap=cm.PuOr, vmin=-maxAct, vmax=maxAct)
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
    
        try:
            tocolor.append((bplot['caps'][i-1], color))
            tocolor.append((bplot['caps'][i], color))
        except:
            pass
        
        try:
            tocolor.append((bplot['fliers'][i-1], color))
            tocolor.append((bplot['fliers'][i], color))
        except:
            pass
        
        try:
            tocolor.append((bplot['whiskers'][i-1], color))
            tocolor.append((bplot['whiskers'][i], color))
        except:
            pass
    
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
    
def plotMapBars(lOrg, title, fname, svg=False, labels=[]):
    '''
    Plot histograms for Kegg mapping statistics
    '''
    plt.clf()
    space = np.array([0.0, 0.2, 0.4, 0.6])
    maxprots = max([x[1] for x in lOrg])
    
    for data in lOrg:
        index = float(lOrg.index(data))
        patch = plt.bar(space + index, np.array(data[1:]),
                width=0.2,
                color=['#3366CC', 'orange', '#D32626', '#cd6464'])
        if index == 0 and len(labels) == 0:
            patch[0].set_label('Size')
            patch[1].set_label('Mapped to Kegg')
            patch[2].set_label('Kegg reactions')
            patch[3].set_label('Unique Kegg reaction IDs')
        elif index == 0 and len(labels) >= 4:
            patch[0].set_label(labels[0])
            patch[1].set_label(labels[1])
            patch[2].set_label(labels[2])
            patch[3].set_label(labels[3])
    
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
    
def plotPanGenomeReactions(conserved, variable, svg=False):
    plt.clf()
    colors=('#D32626','#3366CC')
    patches = plt.pie([conserved, variable], colors=colors,
                      explode=(0.1,0.01),
                      autopct='%1.1f%%',
                      shadow=True)
    plt.legend((patches[0][0], patches[0][1]),
               ('Conserved','Variable'),
               loc=(0,-.1),prop={'size':6})
    plt.title('PanGenome shape (distinct reaction IDs)')
    
    if svg:
        fname = 'pangenome_reaction_shape.svg'
    else:
        fname = 'pangenome_reaction_shape.png'
    
    plt.savefig(fname)
    
    logger.info('PanGenome reactions shape graph saved (%s)'%fname)
    
def isPhenome(project):
    '''
    Checks the presence of phenomic data inside the project
    '''
    return not Biolog(project).isEmpty() 

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
