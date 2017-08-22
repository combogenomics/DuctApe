#!/usr/bin/env python
"""
Database

Storage library

SQLite Database wrappers
"""
# TODO: decorator to catch SQLite exceptions

from ductape.storage.SQLite.dbstrings import dbcreate, dbboost
from ductape.common.utils import get_span
import logging
import sqlite3
import time

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('ductape.database')

################################################################################
# Classes

class Row(object):
    '''
    Class Row
    Holds all the columns as attributes
    Just provide the single row and its description
    '''
    def __init__(self, data, description):
        for field in description:
            try:
                self.__setattr__(field[0],data[description.index(field)])
            except:
                self.__setattr__(field[0],None)

class DBBase(object):
    '''
    Class DB
    General DB handler 
    '''
    def __init__(self, dbname='storage'):
        self.dbname = dbname
        self.connection = None
        self.cursor = None
        self.connect()
    
    def connect(self):
        self.connection = sqlite3.connect(self.dbname)
        
    def getCursor(self):
        if not self.connection:
            self.connect()
        if not self.cursor:
            self.cursor = self.connection.cursor()
    
    def close(self):
        if self.cursor:
            self.cursor.close()
        self.cursor = None
        if self.connection:
            self.connection.close()
        self.connection = None
        
    def create(self):
        '''
        DB creation
        Returns True/False
        '''
        try:
            self.boost()
            
            with self.connection:
                for command in dbcreate.split(';'):
                    self.connection.execute(command+';')
                    
            # Import Biolog data
            b = Biolog(self.dbname)
            b.create()
        except sqlite3.Error as e:
            logger.error('Could not create the database!')
            logger.error(e)
            return False

        return True
    
    def boost(self):
        '''
        The current connection is boosted
        '''
        with self.connection as conn:
            conn.execute(dbboost)
            
    def query(self, sql):
        '''
        Launch a query and returns a generator with each row 
        '''
        with self.connection as conn:
            cursor=conn.execute(sql)
            
        for res in cursor:
            yield Row(res, cursor.description)
    
class Project(DBBase):
    '''
    Class Project
    Handles projects data
    '''
    def __init__(self, dbname='storage'):
        DBBase.__init__(self, dbname)
        
        self.name = None
        self.description = None
        self.kind = None
        self.tmp = None
        self.creation = None
        self.last = None
        self.genome = None
        self.phenome = None
        self.pangenome = None
        
        # Populate the project immediately
        self.getProject()
    
    def __str__(self):
        self.getProject()
        return ' - '.join([
                 str(self.name),
                 str(self.description),
                 str(self.kind),
                 str(self.creation),
                 str(self.last)
                          ])
    
    def isProject(self):
        '''
        Do we have already a project in there?
        '''
        self.getProject()
        if self.name:
            return True
        else:
            return False
    
    def getProject(self):
        '''
        Grep the project informations from the DB
        '''
        # Get the first row (the only one that makes sense)
        with self.connection as conn:
            cursor = conn.execute('select * from project limit 1;')
        
        data = cursor.fetchall()
        if len(data) == 0:
            return
        for field in cursor.description:
            self.__setattr__(field[0],data[0][cursor.description.index(field)])
            
    def addProject(self, name='Project',
                         description='Generic project',
                         kind='generic', tmp=None):
        '''
        If there is no project it is added to db,
        otherwise an exception is thrown
        '''
        if self.isProject():
            logger.warning('Tried to add a project when one is already defined')
            raise Exception('Only one project at a time!')
        
        creation = time.asctime()
        last = time.asctime()
        
        with self.connection as conn:
            conn.execute('''insert into project (`name`, `description`, `kind`,
                                                `tmp`, `creation`, `last`)
                            values (?,?,?,?,?,?);''',
                         (name, description, kind, tmp, creation, last,))
        
    def updateLast(self):
        '''
        Update the last touch timestamp
        '''
        self.getProject()
        with self.connection as conn:
            conn.execute('update project set last = ? where name = ?;',
                         [time.asctime(),self.name,])
        self.getProject()
    
    def setName(self,newName):
        self.getProject()
        with self.connection as conn:
            conn.execute('update project set name = ? where name = ?;',
                         [newName,self.name,])
        # Update the project
        self.getProject()
        
    def setKind(self,newKind):
        self.getProject()
        with self.connection as conn:
            conn.execute('update project set kind = ? where name = ?;',
                         [newKind,self.name,])
        # Update the project
        self.getProject()
    
    def setGenome(self, status):
        '''
        Set the genomic status of the project
        '''
        self.getProject()
        with self.connection as conn:
            conn.execute('update project set genome = ? where name = ?;',
                         [status,self.name,])
        # Update the project
        self.getProject()
    
    def setPhenome(self, status):
        '''
        Set the phenomic status of the project
        '''
        self.getProject()
        with self.connection as conn:
            conn.execute('update project set phenome = ? where name = ?;',
                         [status,self.name,])
        # Update the project
        self.getProject()
        
    def donePanGenome(self):
        '''
        Update the project informing that the pangenome has been done
        '''
        self.getProject()
        with self.connection as conn:
            conn.execute('update project set pangenome = 1 where name = ?;',
                         [self.name,])
        # Update the project
        self.getProject()
        
    def clearPanGenome(self):
        '''
        Update the project informing that the pangenome has been cleared
        '''
        self.getProject()
        with self.connection as conn:
            conn.execute('update project set pangenome = 0 where name = ?;',
                         [self.name,])
        # Update the project
        self.getProject()
        
    def isPanGenome(self):
        '''
        Get the pangenome status
        '''
        self.getProject()
        if self.pangenome == 1:
            return True
        else:
            return False
        
    def setKegg(self, keggver):
        '''
        Set the KEGG database version
        '''
        self.getProject()
        with self.connection as conn:
            conn.execute('update project set kegg = ? where name = ?;',
                         [keggver, self.name,])
        # Update the project
        self.getProject()
        
    def isKegg(self):
        '''
        Check if a Kegg database version has been set
        '''
        self.getProject()
        if not self.kegg:
            return False
        else:
            return True
    
class Organism(DBBase):
    '''
    Class Organism
    Handles the addition and updates on organisms used by the program
    ''' 
    def __init__(self, dbname='storage'):
        DBBase.__init__(self, dbname)
        
    def __len__(self):
        return self.howMany()
    
    def resetProject(self):
        '''
        Reset the project statuses
        '''
        # Reset the project statuses
        oProj = Project(self.dbname)
        oProj.clearPanGenome()
        oProj.setGenome('none')
        oProj.setPhenome('none')
    
    def isOrg(self, org_id):
        '''
        Is this organism already present?
        '''
        with self.connection as conn:
            cursor=conn.execute('select count(*) from organism where org_id=?;',
                                [org_id,])
        return bool(cursor.fetchall()[0][0])
    
    def isMutant(self, org_id):
        '''
        Is this organism a mutant?
        '''
        with self.connection as conn:
            cursor=conn.execute('select mutant from organism where org_id=?;',
                                [org_id,])
        return bool(cursor.fetchall()[0][0])
    
    def getMutReference(self, org_id):
        '''
        Get the reference genome for this mutant
        '''
        query = '''
                select reference
                from organism
                where org_id = ?
                and mutant = 1;
                '''
        
        with self.connection as conn:
            cursor = conn.execute(query, (org_id,))
        
        return cursor.fetchall()[0][0]
    
    def getOrgMutants(self, org_id):
        '''
        Get the list of org_id that are mutants of this organism 
        '''
        query = '''
                select distinct org_id
                from organism
                where reference = ?
                and mutant = 1;
                '''
        
        with self.connection as conn:
            cursor = conn.execute(query, (org_id,))
        
        for mut in cursor:
            yield mut[0]
            
    def howManyMutants(self):
        '''
        Get the overall number of mutants
        '''
        query = '''
                select count(distinct org_id)
                from organism
                where mutant = 1;
                '''
        
        with self.connection as conn:
            cursor = conn.execute(query)
        
        return int(cursor.fetchall()[0][0])
    
    def howMany(self):
        '''
        How many organisms are there?
        '''
        with self.connection as conn:
            cursor=conn.execute('select count(*) from organism;')
        return int(cursor.fetchall()[0][0])
    
    def getAll(self):
        '''
        Returns a list of Row objects about all the organisms
        '''
        with self.connection as conn:
            cursor=conn.execute('select * from organism order by org_id')
            
        for res in cursor:
            yield Row(res, cursor.description)
        
    def getOrg(self, org_id):
        '''
        Get details about one organism
        '''
        if not self.isOrg(org_id):
            return None
        
        with self.connection as conn:
            cursor=conn.execute('select * from organism where org_id=?;',
                                [org_id,])
            
        return Row(cursor.fetchall()[0], cursor.description)
    
    def addOrg(self, org_id, name=None, description=None,
                    orgfile=None, mutant=False, reference=None, mkind='', color=None):
        '''
        Adds a new organism to the db
        Performs some checks on the fields mutant and reference
        If it is a mutant, the reference organism can be null, otherwise
        an exception is raised if it's not present
        '''
        already = self.isOrg(org_id)
        
        mutant = int(mutant)
        if mutant:
            if reference:
                if not self.isOrg(reference):
                    logger.warning('Reference %s is not present yet!'%reference)
                    raise Exception('This reference (%s) is not present yet!'%reference)
        
        with self.connection as conn:
            if not already:
                conn.execute('''insert into organism (`org_id`, `name`,
                                    `description`, `file`, `mutant`, `reference`,
                                    mkind)
                                values (?,?,?,?,?,?,?);''',
                         (org_id, name, description, orgfile, mutant, reference,
                          mkind))
            else:
                conn.execute('''update organism set name = ?,
                                description = ?, mutant = ?, reference = ?,
                                mkind = ? where org_id = ?;''',
                         (name, description, mutant, reference,
                          mkind, org_id))
            
            if color != None:
                conn.execute('''update organism set color = ? where org_id = ?;''',
                         (color, org_id))
        
        if not already:
            # Reset the genomic/phenomic status
            self.setGenomeStatus(org_id, 'none')
            self.setPhenomeStatus(org_id, 'none')
            self.resetProject()
    
    def delAllOrgs(self, cascade=True):
        '''
        Delete all the organsim from the db
        (including all dependent tables)
        '''
        for org in self.getAll():
            self.delOrg(org.org_id, cascade)
            
        self.resetProject()
    
    def delOrg(self, org_id, cascade=False):
        '''
        Delete an organism from the db
        If it's a reference to some mutant remove the children only if cascade 
        '''
        if not self.isOrg(org_id):
            return
        
        with self.connection as conn:
            conn.execute('delete from organism where org_id=?;', (org_id,))
        
        oDel = Genome(self.dbname)
        oBDel = Biolog(self.dbname)
        
        if self.howManyMutants() > 0 and cascade:
            for mut_id in self.getOrgMutants(org_id):
                self.delOrg(mut_id, cascade=True)
                oDel.delProteome(mut_id)
                oBDel.delOrg(mut_id)
                
        oDel.delProteome(org_id)
        oBDel.delOrg(org_id)
        
        self.resetProject()
        
    def setName(self, org_id, newName):
        '''
        Change the name of an organism
        '''
        with self.connection as conn:
            conn.execute('update project set name = ? where org_id = ?;',
                         [newName,org_id,])
            
    def setDescription(self, org_id, newDescr):
        '''
        Change the description of an organism
        '''
        with self.connection as conn:
            conn.execute('update project set description = ? where org_id = ?;',
                         [newDescr,org_id,])
    
    def setFile(self, org_id, newFile):
        '''
        Change the file of an organism
        '''
        with self.connection as conn:
            conn.execute('update project set orgfile = ? where org_id = ?;',
                         [newFile,org_id,])
            
    def setMutant(self, org_id, mutant=1, reference=None):
        '''
        Make this organism no longer a mutant / a mutant: reference can be null
        otherwise an exception is raised if it's not present
        '''
        mutant=int(mutant)
        if reference and mutant:
            if not self.isOrg(reference):
                logger.warning('Reference %s is not present yet!'%reference)
                raise Exception('This reference (%s) is not present yet!'%reference)
        
        with self.connection as conn:
            conn.execute('update organism set mutant = ? where org_id = ?;',
                         [mutant,org_id,])
                         
    def setColor(self, org_id, color):
        '''
        Set the color of this organism (i.e. used for phenomic plots)
        '''
        # Check how many organisms have the same color
        with self.connection as conn:
            cursor=conn.execute('select count(*) from organism where color = ?;',
                        [color,])
        howmany = int(cursor.fetchall()[0][0])
        # We issue just a warning
        if howmany != 0:
            logger.warning('%d organism(s) already use this color (%s)',(howmany, color))
            
        with self.connection as conn:
            conn.execute('update organism set color = ? where org_id = ?;',
                         [color,org_id,])
    
    def resetGenomes(self):
        '''
        Reset each organism genomic status
        '''
        for org in self.getAll():
            self.setGenomeStatus(org.org_id, 'none')
    
    def resetPhenomes(self):
        '''
        Reset each organism phenomic status
        '''
        for org in self.getAll():
            self.setPhenomeStatus(org.org_id, 'none')
    
    def setAllGenomeStatus(self, status):
        '''
        Change the genomic status of all organisms
        '''
        with self.connection as conn:
            conn.execute('update organism set genome = ?;',
                         [status])
    
    def setGenomeStatus(self, org_id, status):
        '''
        Change the genomic status of an organism
        '''
        with self.connection as conn:
            conn.execute('update organism set genome = ? where org_id = ?;',
                         [status,org_id,])
    
    def setAllPhenomeStatus(self, status):
        '''
        Change the genomic status of all organisms
        '''
        with self.connection as conn:
            conn.execute('update organism set phenome = ?;',
                         [status])
    
    def setPhenomeStatus(self, org_id, status):
        '''
        Change the phenomic status of an organism
        '''
        with self.connection as conn:
            conn.execute('update organism set phenome = ? where org_id = ?;',
                         [status,org_id,])
            
class Genome(DBBase):
    '''
    Class Genome
    Handles the addition and updates on Genomic data used by the program
    ''' 
    def __init__(self, dbname='storage'):
        DBBase.__init__(self, dbname)
    
    def resetProject(self):
        '''
        Reset the project Genomic statuses
        '''
        # Reset the project statuses
        oProj = Project(self.dbname)
        oProj.clearPanGenome()
        oProj.setGenome('none')
    
    def updateStatus(self, org_id, status):
        '''
        Update the organism status
        '''
        oOrg = Organism(self.dbname)
        oOrg.setGenomeStatus(org_id, status)
    
    def clearAllGenome(self):
        '''
        Truncate all the tables about the genomic data
        '''
        logger.debug('Clearing genomic data')
        
        with self.connection as conn:
            conn.execute('delete from protein;')
            conn.execute('delete from ortholog;')
            conn.execute('delete from mapko;')
            
        oOrg = Organism(self.dbname)
        oOrg.resetGenomes()
        
        self.resetProject()
    
    def isProt(self, prot_id):
        '''
        Is this protein already present?
        '''
        with self.connection as conn:
            cursor=conn.execute('select count(*) from protein where prot_id=?;',
                                [prot_id,])
        return bool(cursor.fetchall()[0][0])
    
    def areProts(self, prots):
        '''
        Returns False if at least one prot_id is absent
        '''
        with self.connection as conn:
            for prot_id in prots:
                cursor=conn.execute('select count(*) from protein where prot_id=?;',
                                [prot_id,])
                if not bool(cursor.fetchall()[0][0]):
                    logger.warning('Protein %s is not present yet!'%prot_id)
                    return False
        
        return True
    
    def addProteome(self, org_id, pfile):
        '''
        Add a bunch of proteins belonging to org_id (which must be present!)
        The proteins are present in a fasta file, if a particular protein had 
        already been added, no warnings are thrown
        An exception is raised if the org_id is not present in the database
        '''
        from Bio import SeqIO
        
        # Is the organism present?
        oCheck = Organism(self.dbname)
        if not oCheck.isOrg(org_id):
            logger.warning('Organism %s is not present yet!'%org_id)
            raise Exception('This organism (%s) is not present yet!'%org_id)
        
        self.boost()
        
        i = 0
        with self.connection as conn:
            for s in SeqIO.parse(open(pfile),'fasta'):
                conn.execute('insert or replace into protein values (?,?,?,?);',
                         [s.id,org_id,s.description,str(s.seq),])
                i += 1
        
        logger.debug('Added %d protein to organism %s'%(i,org_id))
        
        self.updateStatus(org_id, 'none')
        oProj = Project(self.dbname)
        oProj.clearPanGenome()
             
    def getProt(self, prot_id):
        '''
        Get a specific protein matching the provided prot_id
        '''
        with self.connection as conn:
            cursor=conn.execute('select * from protein where prot_id = ?;',
                                [prot_id,])
        
        data = cursor.fetchall()
        if len(data) == 0:
            return Row([], cursor.description)
        else:
            return Row(data[0], cursor.description)
        
    def getAllProt(self, org_id):
        '''
        Get all the proteins from a specific organism
        '''
        # Is the organism present?
        oCheck = Organism(self.dbname)
        if not oCheck.isOrg(org_id):
            return
        
        with self.connection as conn:
            cursor=conn.execute('select * from protein where org_id = ?;',
                                [org_id,])
            
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getRecords(self, org_id):
        '''
        Get all the proteins from a specific organism
        As SeqRecords objects
        '''
        from Bio import Alphabet
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        
        for prot in self.getAllProt(org_id):
            yield SeqRecord(Seq(prot.sequence,
                                Alphabet.IUPAC.ExtendedIUPACProtein()),
                            id = prot.prot_id)
    
    def howMany(self, org_id):
        '''
        How many proteins for my organism?
        '''
        with self.connection as conn:
            cursor=conn.execute('''select count(*) from protein
                                    where org_id=?;''',[org_id,])
        return int(cursor.fetchall()[0][0])
    
    def delProteome(self, org_id):
        '''
        Remove all the proteins related to a specific organism
        Also the pangenome and the map2ko will be deleted!
        '''
        prots = [x.prot_id for x in self.getAllProt(org_id)]
        self.delKOs(prots)
        
        with self.connection as conn:
            conn.execute('delete from protein where org_id=?;', (org_id,))
            
        self.delPanGenome()
        
        self.resetProject()
        oProj = Project(self.dbname)
        oProj.clearPanGenome()
            
    def addKOs(self, kos, merged=False):
        '''
        Add a bunch of KO IDs mappings
        If merged, the annotation has been taken from the orthology
        '''
        oCheck = Kegg(self.dbname)
        
        self.boost()
        
        for prot_id,ko_id in kos:
            if not self.isProt(prot_id):
                logger.warning('Protein %s is not present yet!'%prot_id)
                raise Exception('This Protein (%s) is not present yet!'%prot_id)
            
            if ko_id.startswith('ko:'):
                ko_id = ko_id.lstrip('ko:')
                
            if not oCheck.isKO('ko:'+ko_id):
                logger.warning('KO %s is not present yet!'%('ko:'+ko_id))
                raise Exception('This KO (%s) is not present yet!'%('ko:'+ko_id))
        
        with self.connection as conn:
            for prot_id,ko_id in kos:
                if ko_id.startswith('ko:'):
                    ko_id = ko_id.lstrip('ko:')
                    
                conn.execute('insert or replace into mapko values (?,?,?);',
                             [prot_id,'ko:'+ko_id,int(merged),])
    
    def getKO(self, prot_id):
        with self.connection as conn:
            cursor=conn.execute('select * from mapko where prot_id = ?;',
                                [prot_id,])
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def delKOs(self, prots):
        with self.connection as conn:
            for prot_id in prots:
                conn.execute('delete from mapko where prot_id=?;', (prot_id,))
                
        self.resetProject()
        
    def delMergedKOs(self):
        with self.connection as conn:
            conn.execute('delete from mapko where indirect=1;')
            
    def howManyMergedKOs(self):
        '''
        Get the number of merged KO links
        '''
        query = '''
                select count(distinct ko_id)
                from mapko m
                where indirect=1;
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query)
        return int(cursor.fetchall()[0][0])
                    
    def addPanGenome(self, orthologs):
        '''
        Add the provided pangenome to the db
        A check on each protein is performed
        An exception is raised if at least one protein is missing
        '''
        # Check if all the proteins are present
        prots = [prot_id for ps in orthologs.values() for prot_id in ps]
        for prot_id in prots:
            if not self.isProt(prot_id):
                logger.warning('Protein %s is not present yet!'%prot_id)
                raise Exception('This Protein (%s) is not present yet!'%prot_id)
        
        self.boost()
        
        # Go for it!
        i = 0
        with self.connection as conn:
            for group_id in orthologs:
                for prot_id in orthologs[group_id]:
                    conn.execute('insert or replace into ortholog values (?,?);',
                             [group_id,prot_id,])
                i += 1
        
        oProj = Project(self.dbname)
        oProj.donePanGenome()
        
        logger.debug('Added %d orthologous groups'%(i))
    
    def getPanGenome(self):
        '''
        Returns a dictionary group_id --> [prot_id, ...]
        '''
        pangenome = {}
        with self.connection as conn:
            cursor = conn.execute('''select * from ortholog;''')
        
        for res in cursor:
            obj = Row(res, cursor.description)
            if obj.group_id not in pangenome:
                pangenome[obj.group_id] = []
            pangenome[obj.group_id].append(obj.prot_id)
            
        return pangenome
        
    def getPanGenomeOrgs(self):
        '''
        Returns a dictionary group_id --> [org_id, ...]
        '''
        pangenome = {}
        with self.connection as conn:
            cursor = conn.execute('''select distinct o.group_id, org_id
                                    from ortholog o, protein p
                                    where o.prot_id = p.prot_id;''')
        
        for res in cursor:
            obj = Row(res, cursor.description)
            if obj.group_id not in pangenome:
                pangenome[obj.group_id] = []
            pangenome[obj.group_id].append(obj.org_id)
        
        for group in pangenome:
            pangenome[group] = sorted(pangenome[group])
            
        return pangenome
    
    def getPanGenomeKOs(self):
        '''
        Returns a double dictionary with pangenome and KO annotations
        group_id --> [prot_id --> [ko_id, ...], ...]
        '''
        ko = {}
        oKegg = Kegg(self.dbname)
        for prot_id, ko_id in oKegg.getAllKO():
            ko[prot_id] = ko.get(prot_id, set())
            ko[prot_id].add(ko_id)
        
        pangenome = self.getPanGenome()
        
        panko = {}
        for group_id, prots in pangenome.iteritems():
            panko[group_id] = {}
            for prot_id in prots:
                panko[group_id][prot_id] = ko.get(prot_id, None)
                
        return panko
        
    def alterPanGenome(self):
        # TODO
        raise NotImplementedError
    
    def _getCore(self):
        '''
        Base method to get the core genome
        '''
        # How many organisms are present?
        oCheck = Organism(self.dbname)
        nOrgs = oCheck.howMany()
        
        query = '''
                select distinct group_id, count(distinct org_id) orgs
                from ortholog o, protein r
                where o.prot_id = r.prot_id
                group by group_id
                HAVING orgs = ?;
                '''
        
        with self.connection as conn:
            cursor = conn.execute(query,
                             [nOrgs,])
        
        return cursor
    
    def getCore(self):
        '''
        Returns a list of orthologous groups names belonging to the Core genome
        '''
        cursor = self._getCore()
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getLenCore(self):
        '''
        Get core genome size
        '''
        cursor = self._getCore()
        
        i = 0
        for res in cursor:
            i += 1
            
        return i
    
    def _getDisp(self):
        '''
        Base method to get the dispensable genome
        '''
        # How many organisms are present?
        oCheck = Organism(self.dbname)
        nOrgs = oCheck.howMany()
        
        query = '''
                select distinct group_id, count(distinct org_id) orgs
                from ortholog o, protein r
                where o.prot_id = r.prot_id
                group by group_id
                HAVING orgs < ?;
                '''
        
        with self.connection as conn:
            cursor = conn.execute(query,
                             [nOrgs,])
        
        return cursor
    
    def getDisp(self):
        '''
        Returns a list of orthologous groups names belonging to the Dispensable genome
        '''
        cursor = self._getDisp()
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getLenDisp(self):
        '''
        Get dispensable genome size
        '''
        cursor = self._getDisp()
        
        i = 0
        for res in cursor:
            i += 1
            
        return i
    
    def _getAcc(self):
        '''
        Base method to get the accessory genome
        '''
        # How many organisms are present?
        oCheck = Organism(self.dbname)
        nOrgs = oCheck.howMany()
        
        query = '''
                select distinct group_id, count(distinct org_id) orgs
                from ortholog o, protein r
                where o.prot_id = r.prot_id
                group by group_id
                HAVING orgs < ?
                and orgs > ?;
                '''
        
        with self.connection as conn:
            cursor = conn.execute(query,
                             [nOrgs,1,])
            
        return cursor
    
    def getAcc(self):
        '''
        Returns a list of orthologous groups names belonging to the Accessory genome
        '''
        cursor = self._getAcc()
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getLenAcc(self):
        '''
        Get accessory genome size
        '''
        cursor = self._getAcc()
        
        i = 0
        for res in cursor:
            i += 1
            
        return i
    
    def _getUni(self):
        '''
        Base method to get the unique genome
        '''
        query = '''
                select distinct group_id, count(distinct org_id) orgs
                from ortholog o, protein r
                where o.prot_id = r.prot_id
                group by group_id
                HAVING orgs = ?;
                '''
        
        with self.connection as conn:
            cursor = conn.execute(query,
                             [1,])
        
        return cursor
    
    def getUni(self):
        '''
        Returns a list of orthologous groups names belonging to the Unique genome
        '''
        cursor = self._getUni()
        
        for res in cursor:
            yield Row(res, cursor.description)    
            
    def getLenUni(self):
        '''
        Get unique genome size
        '''
        cursor = self._getUni()
        
        i = 0
        for res in cursor:
            i += 1
            
        return i  
    
    def getGroupNum(self, group_id):
        '''
        Get the number of organisms having the provided ortholog
        '''
        query = '''
                select count(distinct org_id) num
                from ortholog o, protein p
                where o.prot_id = p.prot_id
                and group_id = ?;
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query,[group_id,])
        
        return int(Row(cursor.fetchall()[0], cursor.description).num)
        
    def delPanGenome(self):
        '''
        Remove all the entries about the pangenome
        '''
        with self.connection as conn:
            conn.execute('delete from ortholog')
            
        self.resetProject()
        oProj = Project(self.dbname)
        oProj.clearPanGenome()

class Kegg(DBBase):
    '''
    Class Kegg
    Handles all the data about Kegg entries
    '''
    def __init__(self, dbname='storage'):
        DBBase.__init__(self, dbname)
    
    def clear(self):
        '''
        Remove all the KEGG data
        '''
        self.boost()
        
        # Delete the data
        with self.connection as conn:
            conn.execute('delete from ko;')
            conn.execute('delete from ko_react;')
            conn.execute('delete from compound;')
            conn.execute('delete from pathway;')
            conn.execute('delete from reaction;')
            conn.execute('delete from rpair;')
            conn.execute('delete from comp_path;')
            conn.execute('delete from react_comp;')
            conn.execute('delete from react_path;')
            conn.execute('delete from rpair_react;')
        
        # "Update" the release number
        proj = Project(self.dbname)
        proj.setKegg(None)
    
    def exportKegg(self):
        '''
        Generator for kegg data export
        All the relevant data for the kegg db is extracted
        '''
        # Release?
        try:
            oCheck = Project(self.dbname)
            if oCheck.isKegg():
                yield '\t'.join(['release', str(oCheck.kegg)])
            else:
                yield '\t'.join(['release', str(None)])
        except:
            # Testing bugfix for old DBs
            yield '\t'.join(['release', str(None)])
        
        with self.connection as conn:
            cursor=conn.execute('select * from ko;')
            
        for res in cursor:
            yield '\t'.join(['ko'] + 
                            [str(res[cursor.description.index(field)]) 
                             for field in cursor.description])
        
        with self.connection as conn:
            cursor=conn.execute('select * from compound;')
            
        for res in cursor:
            yield '\t'.join(['compound'] + 
                            [str(res[cursor.description.index(field)]) 
                             for field in cursor.description])
    
        with self.connection as conn:
            conn.text_factory = str
            
            cursor=conn.execute('select * from pathway;')
            
        # Exceptionally ugly exception
        # Newlines chars in html maps are converted in DUCTAPENEWLINEHERE
        # Which is stupid, but for know it will suffice
        for res in cursor:
            yield '\t'.join(['pathway'] + 
                            [str(res[cursor.description.index(field)]).replace('\n','DUCTAPENEWLINEHERE').replace('\t','  ') 
                             for field in cursor.description])
    
        with self.connection as conn:
            cursor=conn.execute('select * from reaction;')
            
        for res in cursor:
            yield '\t'.join(['reaction'] + 
                            [str(res[cursor.description.index(field)]) 
                             for field in cursor.description])
            
        with self.connection as conn:
            cursor=conn.execute('select * from rpair;')
            
        for res in cursor:
            yield '\t'.join(['rpair'] + 
                            [str(res[cursor.description.index(field)]) 
                             for field in cursor.description])
        
        with self.connection as conn:
            cursor=conn.execute('select * from ko_react;')
            
        for res in cursor:
            yield '\t'.join(['ko_react'] + 
                            [str(res[cursor.description.index(field)]) 
                             for field in cursor.description])
        
        with self.connection as conn:
            cursor=conn.execute('select * from comp_path;')
            
        for res in cursor:
            yield '\t'.join(['comp_path'] + 
                            [str(res[cursor.description.index(field)]) 
                             for field in cursor.description])
            
        with self.connection as conn:
            cursor=conn.execute('select * from react_comp;')
            
        for res in cursor:
            yield '\t'.join(['react_comp'] + 
                            [str(res[cursor.description.index(field)]) 
                             for field in cursor.description])
    
        with self.connection as conn:
            cursor=conn.execute('select * from react_path;')
            
        for res in cursor:
            yield '\t'.join(['react_path'] + 
                            [str(res[cursor.description.index(field)]) 
                             for field in cursor.description])
            
        with self.connection as conn:
            cursor=conn.execute('select * from rpair_react;')
            
        for res in cursor:
            yield '\t'.join(['rpair_react'] + 
                            [str(res[cursor.description.index(field)]) 
                             for field in cursor.description])
    
    def importKegg(self, infile):
        '''
        Imports the content of the file object inside the kegg tables
        In case of errors there should be a rollback
        '''
        self.boost()
        
        with self.connection as conn:
            conn.text_factory = str
            
            for l in infile:
                if l.lstrip().startswith('#'):continue
            
                s = l.rstrip('\n').split('\t')
                # Special case #1
                if s[0] == 'release':
                    if s[1] == 'None':
                        release = None
                    else:
                        release = s[1]
                else:
                    # Special case #2
                    if s[0] == 'pathway':
                        for i in range(len(s)):
                            s[i] = s[i].replace('DUCTAPENEWLINEHERE','\n')
                    
                    for i in range(len(s)):
                        if s[i] == 'None':
                            s[i] = None
                    
                    values = ''
                    for i in range(len(s[1:])):
                        values += '''?, '''
                    values = values.rstrip(', ')
                    query = '''insert or replace into %s values (%s);'''%(s[0], values)
                    
                    conn.execute(query, [str(x) for x in s[1:]])
        
        # Last step
        if release:
            proj = Project(self.dbname)
            proj.setKegg(release)            
    
    def addDraftKOs(self, ko):
        '''
        Add new KOs (skipping if they are already present)
        the input is a list, so no details about this KOs are there yet
        '''
        self.boost()
        
        with self.connection as conn:
            for ko_id in ko:
                conn.execute('insert or ignore into ko (`ko_id`) values (?);',
                     ('ko:'+ko_id,))
    
    def addKOs(self, ko):
        '''
        Add new KOs (ignoring errors if they are already present)
        the input is a dictionary
        ko_id --> name, description
        '''
        self.boost()
        
        with self.connection as conn:
            for ko_id,values in ko.iteritems():
                name = values[0]
                if len(values) > 1:
                    description = values[1]
                else:
                    description = ''
                conn.execute('insert or replace into ko values (?,?,?,?);',
                     (ko_id,name,description,1,))
    
    def addKOReacts(self, koreact):
        '''
        An exception is thrown if such IDs are not present
        '''
        for ko_id in koreact:
            if not self.isKO(ko_id):
                logger.warning('KO %s is not present yet!'
                                %ko_id)
                raise Exception('This KO (%s) is not present yet!'
                                %ko_id)
            for re_id in koreact[ko_id]:
                if not self.isReaction(re_id):
                    logger.warning('Reaction %s is not present yet!'
                                %re_id)
                    raise Exception('This reaction (%s) is not present yet!'
                                %re_id)
        
        self.boost()
        
        with self.connection as conn:
            for ko_id in koreact:
                for re_id in koreact[ko_id]:
                    conn.execute('insert or ignore into ko_react values (?,?);',
                                 (ko_id,re_id,))
    
    def getKO2Analyze(self):
        '''
        Get all the ko_id to be analyzed
        '''
        with self.connection as conn:
            cursor=conn.execute('select ko_id from ko where analyzed = 0;')
            
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getAllIDs(self):
        '''
        Get all the Kegg IDs already analyzed
        '''
        with self.connection as conn:
            cursor=conn.execute('select ko_id from ko where analyzed = 1;')
            
        for res in cursor:
            yield res[0]
            
        with self.connection as conn:
            cursor=conn.execute('select re_id from reaction;')
            
        for res in cursor:
            yield res[0]
            
        with self.connection as conn:
            cursor=conn.execute('select co_id from compound;')
            
        for res in cursor:
            yield res[0]
    
        with self.connection as conn:
            cursor=conn.execute('select path_id from pathway;')
            
        for res in cursor:
            yield res[0]
            
        with self.connection as conn:
            cursor=conn.execute('select rp_id from rpair;')
            
        for res in cursor:
            yield res[0]
            
    def getKO(self, ko_id):
        '''
        Get a specific ko_id
        '''
        if not self.isKO(ko_id):
            return None
        
        with self.connection as conn:
            cursor=conn.execute('select * from ko where ko_id=?;',
                                (ko_id,))
            
        return Row(cursor.fetchall()[0], cursor.description)
    
    def getAllKO(self, org_id=None, merged=False):
        '''
        Get all the prot_id, ko_id pair iterator from a specific org_id
        If org_id is not provided, all the pairs are provided
        If merged is True, only those proteins annotated by orthology are provided
        '''
        if not org_id:
            if not merged:
                query = '''
                        select distinct prot_id, ko_id
                        from mapko
                        '''
            else:
                query = '''
                        select distinct prot_id, ko_id
                        from mapko
                        where indirect=1
                        '''
                
            with self.connection as conn:
                cursor=conn.execute(query)
        
        else:
            if not merged:
                query = '''
                        select distinct m.prot_id, ko_id
                        from mapko m, protein p
                        where m.prot_id = p.prot_id
                        and org_id = ?;
                        '''
            else:
                query = '''
                        select distinct m.prot_id, ko_id
                        from mapko m, protein p
                        where m.prot_id = p.prot_id
                        and org_id = ?
                        and indirect=1;
                        '''
        
            with self.connection as conn:
                cursor=conn.execute(query,[org_id,])
            
        for res in cursor:
            yield res[0], res[1]
    
    def getMultipleKOs(self, org_id=None):
        '''
        Generator to the multiple annotated proteins
        If org_id is not provided, all the pairs are provided
        '''
        if not org_id:
            query = '''
                    select prot_id, count(*)
                    from mapko
                    group by prot_id
                    having count(*)>1
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query)
            
        else:
            query = '''
                    select m.prot_id, count(*)
                    from mapko m, protein p
                    where m.prot_id=p.prot_id
                    and org_id=?
                    group by m.prot_id
                    having count(*)>1
                    '''
        
            with self.connection as conn:
                cursor=conn.execute(query,[org_id,])
            
        for res in cursor:
            query = '''
                    select distinct ko_id
                    from mapko
                    where prot_id = ?;
                    '''
            
            with self.connection as conn:
                cursor1=conn.execute(query,[res[0]])
                
            for res1 in cursor1:
                yield res[0], res1[0]
    
    def isKO(self, ko_id):
        '''
        Is this ko_id already present?
        '''
        try:
            with self.connection as conn:
                cursor=conn.execute('select count(*) from ko where ko_id=?;',
                                    (ko_id,))
            return bool(cursor.fetchall()[0][0])
        except Exception as e:
            logger.debug('Got error %s on id %s, assuming id is present'%
                         (str(e),ko_id))
            return True
    
    def addReactions(self, react):
        '''
        Add new reactions (ignoring errors if they are already present)
        the input is a dictionary
        re_id --> name, description, enzyme
        '''
        self.boost()
        
        with self.connection as conn:
            for re_id, values in react.iteritems():
                name = values[0]
                
                if len(values) > 1:
                    description = values[1]
                else:
                    description = None
                    
                if len(values) > 2:
                    enzyme = values[2]
                else:
                    enzyme = None
                    
                conn.execute('insert or ignore into reaction values (?,?,?,?);',
                     (re_id,name,description,enzyme,))
    
    def isRPair(self, rp_id):
        '''
        Is this rp_id already present?
        '''
        try:
            with self.connection as conn:
                cursor=conn.execute('select count(*) from rpair where rp_id=?;',
                                    (rp_id,))
            return bool(cursor.fetchall()[0][0])
        except Exception as e:
            logger.debug('Got error %s on id %s, assuming id is present'%
                         (str(e),rp_id))
            return True
        
    def hasRPairMain(self, path_id):
        '''
        Inspect if a pathway as at least one rpair of kind "main"
        '''
        query = '''
                select count(*)
                from rpair rp, rpair_react r,  reaction re, react_path p
                where r.rp_id = rp.rp_id
                and r.re_id = re.re_id
                and r.re_id=p.re_id
                and path_id=?
                and kind like "%main%";
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query,
                                (path_id,))
        return bool(cursor.fetchall()[0][0])
    
    def getAllRPairCoPaths(self):
        '''
        Generator to compounds to pathways association
        Only those compounds partecipating in a "main" RPair reaction are returned
        '''
        with self.connection as conn:
            cursor=conn.execute('''select distinct cp.co_id, cp.path_id
                                    from rpair rp, comp_path cp
                                    where kind like "%main%"
                            and (rp.co1 = cp.co_id or rp.co2 = cp.co_id);''')
            
        for res in cursor:
            yield Row(res, cursor.description)
       
    def getAllRPairsReacts(self, org_id=None, path_id=None):
        '''
        Generator to single reactions with main rpairs
        If org_id is set, the organism specific subset is retrieved.
        If path_id is set, only those reactiomns from the desired pathway will be retrieved
        '''
        with self.connection as conn:
            if not org_id:
                if not path_id:
                    cursor=conn.execute('''select distinct r.re_id, co1, co2, re.name
                                    from rpair rp, rpair_react r,  reaction re
                                    where r.rp_id = rp.rp_id
                                    and r.re_id = re.re_id
                                    and kind like "%main%";''')
                else:
                    cursor=conn.execute('''select distinct r.re_id, co1, co2, re.name
                                    from rpair rp, rpair_react r,  reaction re, react_path p
                                    where r.rp_id = rp.rp_id
                                    and r.re_id = re.re_id
                                    and r.re_id=p.re_id
                                    and path_id=?
                                    and kind like "%main%";''',
                                    [path_id,])
            else:
                if not path_id:
                    cursor=conn.execute('''select distinct k.re_id, co1, co2, re.name
                                    from ko_react k, mapko m, protein p, rpair_react rr, rpair rp, reaction re
                                    where k.ko_id = m.ko_id
                                    and p.prot_id = m.prot_id
                                    and k.re_id = rr.re_id
                                    and rr.rp_id = rp.rp_id
                                    and rr.re_id=re.re_id
                                    and kind like "%main%"
                                    and org_id = ?;''',
                                    [org_id,])
                else:
                    cursor=conn.execute('''select distinct k.re_id, co1, co2, re.name
                                    from ko_react k, mapko m, protein p, rpair_react rr, rpair rp, reaction re, react_path p
                                    where k.ko_id = m.ko_id
                                    and p.prot_id = m.prot_id
                                    and k.re_id = rr.re_id
                                    and rr.rp_id = rp.rp_id
                                    and rr.re_id=re.re_id
                                    and kind like "%main%"
                                    and org_id = ?
                                    and re.re_id=p.re_id
                                    and path_id = ?;''',
                                    [org_id,path_id,])
            
        for res in cursor:
            yield Row(res, cursor.description)
            
    def addRPairs(self, rp):
        '''
        Add new reactions (ignoring errors if they are already present)
        the input is a dictionary
        co_id --> name, description
        '''
        self.boost()
        
        with self.connection as conn:
            for rp_id, values in rp.iteritems():
                co1, co2, kind = values[0], values[1], values[2]
                conn.execute('insert or ignore into rpair values (?,?,?,?);',
                     (rp_id,co1,co2,kind,))
    
    def addReactRPairs(self, reactrpair):
        '''
        An exception is thrown if such IDs are not present
        '''
        for re_id in reactrpair:
            if not self.isReaction(re_id):
                logger.warning('Reaction %s is not present yet!'
                                %re_id)
                raise Exception('This reaction (%s) is not present yet!'
                                %re_id)
            for rp_id in reactrpair[re_id]:
                if not self.isRPair(rp_id):
                    logger.warning('RPair %s is not present yet!'
                                %rp_id)
                    raise Exception('This rpair (%s) is not present yet!'
                                %rp_id)
        
        self.boost()
        
        with self.connection as conn:
            for re_id in reactrpair:
                for rp_id in reactrpair[re_id]:
                    conn.execute('insert or ignore into rpair_react values (?,?);',
                                 (rp_id,re_id,))
                    
    def addRPairReacts(self, rpairreact):
        '''
        An exception is thrown if such IDs are not present
        '''
        for rp_id in rpairreact:
            if not self.isRPair(rp_id):
                logger.warning('RPair %s is not present yet!'
                                %rp_id)
                raise Exception('This rpair (%s) is not present yet!'
                                %rp_id)
            for re_id in rpairreact[rp_id]:
                if not self.isReaction(re_id):
                    logger.warning('Reaction %s is not present yet!'
                                %re_id)
                    raise Exception('This reaction (%s) is not present yet!'
                                %re_id)
        
        self.boost()
        
        with self.connection as conn:
            for rp_id in rpairreact:
                for re_id in rpairreact[rp_id]:
                    conn.execute('insert or ignore into rpair_react values (?,?);',
                                 (rp_id,re_id,))
    
    def addReactComps(self, reactcomp):
        '''
        An exception is thrown if such IDs are not present
        '''
        for re_id in reactcomp:
            if not self.isReaction(re_id):
                logger.warning('Reaction %s is not present yet!'
                                %re_id)
                raise Exception('This reaction (%s) is not present yet!'
                                %re_id)
            for co_id in reactcomp[re_id]:
                if not self.isCompound(co_id):
                    logger.warning('Compound %s is not present yet!'
                                %co_id)
                    raise Exception('This compound (%s) is not present yet!'
                                %co_id)
        
        self.boost()
        
        with self.connection as conn:
            for re_id in reactcomp:
                for co_id in reactcomp[re_id]:
                    conn.execute('insert or ignore into react_comp values (?,?);',
                                 (re_id,co_id,))
                    
    def addCompReacts(self, compreact):
        '''
        An exception is thrown if such IDs are not present
        '''
        for co_id in compreact:
            if not self.isCompound(co_id):
                logger.warning('Compound %s is not present yet!'
                                %co_id)
                raise Exception('This compound (%s) is not present yet!'
                                %co_id)
            for re_id in compreact[co_id]:
                if not self.isReaction(re_id):
                    logger.warning('Reaction %s is not present yet!'
                                %re_id)
                    raise Exception('This reaction (%s) is not present yet!'
                                %re_id)
        
        self.boost()
        
        with self.connection as conn:
            for co_id in compreact:
                for re_id in compreact[co_id]:
                    conn.execute('insert or ignore into react_comp values (?,?);',
                                 (re_id,co_id,))
    
    def getReaction(self, re_id):
        if not self.isReaction(re_id):
            return None
        
        with self.connection as conn:
            cursor=conn.execute('select * from reaction where re_id=?;',
                                [re_id,])
            
        return Row(cursor.fetchall()[0], cursor.description)
    
    def getAllReactions(self, org_id):
        '''
        Get all the prot_id, re_id pair iterator from a specific org_id
        '''
        query = '''
                select distinct m.prot_id, r.re_id, name, r1.description
                from mapko m, protein p, ko_react r, reaction r1
                where m.prot_id = p.prot_id
                and m.ko_id = r.ko_id
		        and r.re_id=r1.re_id
                and org_id = ?;
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query,[org_id,])
            
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getAllECNumbers(self, org_id):
        '''
        Get all the prot_id, EC numbers pair iterator from a specific org_id
        '''
        query = '''
                select distinct m.prot_id, enzyme
                from mapko m, protein p, ko_react r, reaction r1
                where m.prot_id = p.prot_id
                and m.ko_id = r.ko_id
                and r.re_id = r1.re_id
                and org_id = ?
                and enzyme not NULL
                and enzyme != '';
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query,[org_id,])
            
        for res in cursor:
            yield res[0], res[1]
    
    def isReaction(self, re_id):
        '''
        Is this re_id already present?
        '''
        try:
            with self.connection as conn:
                cursor=conn.execute('select count(*) from reaction where re_id=?;',
                                    (re_id,))
            return bool(cursor.fetchall()[0][0])
        except Exception as e:
            logger.debug('Got error %s on id %s, assuming id is present'%
                         (str(e),re_id))
            return True
        
    def addCompounds(self, co):
        '''
        Add new reactions (ignoring errors if they are already present)
        the input is a dictionary
        co_id --> name, description
        '''
        self.boost()
        
        with self.connection as conn:
            for co_id, values in co.iteritems():
                name = values[0]
                if len(values) > 1:
                    description = values[1]
                else:
                    description = ''
                conn.execute('insert or ignore into compound values (?,?,?);',
                     (co_id,name,description,))
    
    def getCompound(self, co_id):
        if not self.isCompound(co_id):
            return None
        
        with self.connection as conn:
            cursor=conn.execute('select * from compound where co_id=?;',
                                [co_id,])
            
        return Row(cursor.fetchall()[0], cursor.description)
    
    def getAllCompounds(self, path_id=None):
        '''
        Generator to all the compounds
        If path_id is set, only those compounds from the desired pathway will be retrieved
        '''
        with self.connection as conn:
            if not path_id:
                cursor=conn.execute('''select * from compound;''')
            else:
                cursor=conn.execute('''
                                    select *
                                    from compound c, comp_path p
                                    where c.co_id = p.co_id
                                    and path_id = ?;
                                    ''',
                                    [path_id,])
            
        for res in cursor:
            yield Row(res, cursor.description)
    
    def isCompound(self, co_id):
        '''
        Is this co_id already present?
        '''
        try:
            with self.connection as conn:
                cursor=conn.execute('select count(*) from compound where co_id=?;',
                                    (co_id,))
            return bool(cursor.fetchall()[0][0])
        except Exception as e:
            logger.debug('Got error %s on id %s, assuming id is present'%
                         (str(e),co_id))
            return True
    
    def addPathways(self, path):
        '''
        Add new pathways (ignoring errors if they are already present)
        the input is a dictionary
        path_id --> name, description
        '''
        self.boost()
        
        with self.connection as conn:
            for path_id, values in path.iteritems():
                name = values[0]
                if len(values) > 1:
                    description = values[1]
                else:
                    description = ''
                conn.execute('''insert or replace into pathway
                                (path_id,name,description) values (?,?,?);''',
                     (path_id,name,description,))
    
    def addPathHtml(self, path):
        '''
        Add pathways HTML (ignoring if the pathways are not present)
        the input is a dictionary
        path_id --> html
        '''
        self.boost()
        
        with self.connection as conn:
            conn.text_factory = str
            
            for path_id, html in path.iteritems():
                if not html:continue
                html = '\n'.join(html)
                conn.execute('''update pathway set html=? where path_id=?;''',
                     (html,path_id,))
                
    def addPathReacts(self, pathreact):
        '''
        An exception is thrown if such IDs are not present
        '''
        for path_id in pathreact:
            if not self.isPathway(path_id):
                logger.warning('Pathway %s is not present yet!'
                                %path_id)
                raise Exception('This pathway (%s) is not present yet!'
                                %path_id)
            for re_id in pathreact[path_id]:
                if not self.isReaction(re_id):
                    logger.warning('Reaction %s is not present yet!'
                                %re_id)
                    raise Exception('This reaction (%s) is not present yet!'
                                %re_id)
        
        self.boost()
        
        with self.connection as conn:
            for path_id in pathreact:
                for re_id in pathreact[path_id]:
                    conn.execute('insert or ignore into react_path values (?,?);',
                                 (re_id,path_id,))
    
    def getReactPath(self, re_id):
        '''
        Get all the pathways linked to a reaction ID
        '''
        query = '''select distinct path_id
                   from react_path
                   where re_id = ?'''
        
        with self.connection as conn:
            cursor=conn.execute(query, [re_id,])
            
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getPathReacts(self):
        '''
        Get all the pathway - reaction links
        '''
        query = '''select path_id, re_id
                   from react_path
                   order by path_id'''
        
        with self.connection as conn:
            cursor=conn.execute(query)
            
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getPathComps(self, path_id=None):
        '''
        Get all the pathway - compounds links
        If path_id is set, only those compounds from the desired path are retrieved
        '''
        with self.connection as conn:
            if not path_id:
                cursor=conn.execute('''select path_id, co_id
                                       from comp_path
                                       order by path_id;''')
            else:
                cursor=conn.execute('''select path_id, co_id
                                       from comp_path
                                       where path_id = ?
                                       order by co_id;''',[path_id,])
            
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getCompPaths(self, co_id):
        '''
        Get the pathways related to the provided co_id
        '''
        with self.connection as conn:
            cursor=conn.execute('''select path_id, co_id
                                   from comp_path
                                   where co_id = ?
                                   order by path_id;''',[co_id,])
            
        for res in cursor:
            yield Row(res, cursor.description)
                    
    def addPathComps(self, pathcomp):
        '''
        An exception is thrown if such IDs are not present
        '''
        for path_id in pathcomp:
            if not self.isPathway(path_id):
                logger.warning('Pathway %s is not present yet!'
                                %path_id)
                raise Exception('This pathway (%s) is not present yet!'
                                %path_id)
            for co_id in pathcomp[path_id]:
                if not self.isCompound(co_id):
                    logger.warning('Compound %s is not present yet!'
                                %co_id)
                    raise Exception('This compound (%s) is not present yet!'
                                %co_id)
        
        self.boost()
        
        with self.connection as conn:
            for path_id in pathcomp:
                for co_id in pathcomp[path_id]:
                    conn.execute('insert or ignore into comp_path values (?,?);',
                                 (co_id,path_id,))
                    
    def addPathMaps(self, pathmap):
        '''
        An exception is thrown if such ID is not present
        '''
        for path_id in pathmap:
            if not self.isPathway(path_id):
                logger.warning('Pathway %s is not present yet!'
                                %path_id)
                raise Exception('This pathway (%s) is not present yet!'
                                %path_id)
        
        self.boost()
        
        with self.connection as conn:
            conn.text_factory = str
            
            for path_id in pathmap:
                conn.execute('insert or ignore into pathmap (path_id, html) values (?,?);',
                                 (path_id,'\n'.join(pathmap[path_id]),))
    
    def addPathPics(self, pathpic):
        '''
        An exception is thrown if such ID is not present
        '''
        for path_id in pathpic:
            if not self.isPathway(path_id):
                logger.warning('Pathway %s is not present yet!'
                                %path_id)
                raise Exception('This pathway (%s) is not present yet!'
                                %path_id)
        
        self.boost()
        
        with self.connection as conn:
            for path_id in pathpic:
                pic = open(pathpic[path_id])
                conn.execute('update pathmap set png = ? where path_id = ?;',
                                 (sqlite3.Binary(pic.read()),path_id,))
    
    def getPathway(self, path_id):
        if not self.isPathway(path_id):
            return None
        
        with self.connection as conn:
            conn.text_factory = str
            
            cursor=conn.execute('select * from pathway where path_id=?;',
                                [path_id,])
            
        return Row(cursor.fetchall()[0], cursor.description)
    
    def getAllPathways(self, onlymain=False):
        '''
        Generator to the single pathways
        If onlymain is set, only pathways with RPairs main are returned
        '''
        
        query = '''select * from pathway order by path_id;'''
        
        with self.connection as conn:
            conn.text_factory = str
            cursor=conn.execute(query)
            
        for res in cursor:
            path = Row(res, cursor.description)
            if not self.hasRPairMain(path.path_id):
                continue
            yield path
    
    def isPathway(self, path_id):
        '''
        Is this path_id already present?
        '''
        try:
            with self.connection as conn:
                cursor=conn.execute('select count(*) from pathway where path_id=?;',
                                    [path_id,])
            return bool(cursor.fetchall()[0][0])
        except Exception as e:
            logger.debug('Got error %s on id %s, assuming id is present'%
                         (str(e),path_id))
            return True
    
    def getReactNum(self, re_id):
        '''
        Get the number of organisms sharing this re_id
        '''
        query = '''
                select count(distinct org_id) num
                from ortholog o, mapko m, protein p, ko_react k
                where o.prot_id = p.prot_id
                and p.prot_id = m.prot_id
                and m.ko_id = k.ko_id
                and re_id=?
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query,[re_id,])
        
        return int(Row(cursor.fetchall()[0], cursor.description).num)
    
    def getAllReactNum(self):
        '''
        Get the number of organisms sharing each re_id
        '''
        query = '''
                select k.re_id, count(distinct org_id) num
                from ortholog o, mapko m, protein p, ko_react k
                where o.prot_id = p.prot_id
                and p.prot_id = m.prot_id
                and m.ko_id = k.ko_id
                group by k.re_id
                order by num
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query)
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getMappedRPairsReact(self, path_id=None):
        '''
        Get all the RPairs Reacts in the pangenome
        '''
        if not path_id:
            query = '''
                    select distinct k.re_id, co1, co2, re.name
                    from ko_react k, mapko m, protein p, rpair_react rr, rpair rp, reaction re
                    where k.ko_id = m.ko_id
                    and p.prot_id = m.prot_id
                    and k.re_id = rr.re_id
                    and rr.rp_id = rp.rp_id
                    and rr.re_id=re.re_id
                    and kind like "%main%"
                    '''
        else:
            query = '''
                    select distinct k.re_id, co1, co2, re.name
                    from ko_react k, mapko m, protein p, rpair_react rr, rpair rp, reaction re, react_path p1
                    where k.ko_id = m.ko_id
                    and p.prot_id = m.prot_id
                    and k.re_id = rr.re_id
                    and rr.rp_id = rp.rp_id
                    and rr.re_id=re.re_id
                    and re.re_id=p1.re_id
                    and p1.path_id=?
                    and kind like "%main%"
                    '''
    
        with self.connection as conn:
            if not path_id:
                cursor=conn.execute(query)
            else:
                cursor=conn.execute(query,
                                    [path_id,])
            
        rnums = {}
        for r in self.getAllReactNum():
            rnums[r.re_id] = r.num
            
        for res in cursor:
            r = Row(res, cursor.description)
            setattr(r, 'weight', rnums[r.re_id])
            yield r
    
    def getExclusiveRPairsReact(self, path_id=None):
        '''
        Get all the exclusive RPairs Reacts in the pangenome
        core, dispensable, accessory, unique
        note: the dispensable genome includes the accessory and the unique
        '''
        genome = Genome(self.dbname)
        
        if not path_id:
            query = '''
                    select distinct k.re_id, co1, co2, re.name, o.group_id
                        from ko_react k, mapko m, rpair_react rr, rpair rp, reaction re, ortholog o
                        where k.ko_id = m.ko_id
                        and m.prot_id = o.prot_id
                        and k.re_id = rr.re_id
                        and rr.rp_id = rp.rp_id
                        and rr.re_id=re.re_id
                        and kind like "%main%"
                    '''
        else:
            query = '''
                    select distinct k.re_id, co1, co2, re.name, o.group_id
                    from ko_react k, mapko m, rpair_react rr, rpair rp, reaction re, ortholog o, react_path p1
                    where k.ko_id = m.ko_id
                    and m.prot_id = o.prot_id
                    and k.re_id = rr.re_id
                    and rr.rp_id = rp.rp_id
                    and rr.re_id=re.re_id
                    and re.re_id=p1.re_id
                    and p1.path_id=?
                    and kind like "%main%"
                    '''
        
        with self.connection as conn:
            if not path_id:
                cursor=conn.execute(query)
            else:
                cursor=conn.execute(query,[path_id,])
        
        ocore = set([x.group_id for x in genome.getCore()])
        odisp = set([x.group_id for x in genome.getDisp()])
        oacc = set([x.group_id for x in genome.getAcc()])
        ouni = set([x.group_id for x in genome.getUni()])
            
        rall = [Row(res, cursor.description) for res in cursor]

        rcore = set([r for r in filter(lambda x: x.group_id in ocore, rall)])
        rdisp = set([r for r in filter(lambda x: x.group_id in odisp, rall)])
        racc = set([r for r in filter(lambda x: x.group_id in oacc, rall)])
        runi = set([r for r in filter(lambda x: x.group_id in ouni, rall)])
        
        rnums = {}
        for r in self.getAllReactNum():
            rnums[r.re_id] = r.num
        
        core = {}
        for r in rcore:
            setattr(r, 'weight', rnums[r.re_id])
            core[r.re_id+r.co1+r.co2] = r
        disp = {}
        for r in rdisp:
            setattr(r, 'weight', rnums[r.re_id])
            disp[r.re_id+r.co1+r.co2] = r
        acc = {}
        for r in racc:
            setattr(r, 'weight', rnums[r.re_id])
            acc[r.re_id+r.co1+r.co2] = r
        uni = {}
        for r in runi:
            setattr(r, 'weight', rnums[r.re_id])
            uni[r.re_id+r.co1+r.co2] = r
        
        ecore = set(core.keys()).difference(set(disp.keys()))
        edisp = set(disp.keys()).difference(set(core.keys()))
        eacc = set(acc.keys()).difference(set(core.keys()), set(uni.keys()))
        euni = set(uni.keys()).difference(set(core.keys()), set(acc.keys()))
        
        return (set([core[x] for x in ecore]),
                set([disp[x] for x in edisp]),
                set([acc[x] for x in eacc]),
                set([uni[x] for x in euni]))
    
    def getExclusiveRPairsReactMutants(self, ref_id, muts=set(), path_id=None):
        '''
        Return the RPairs Reacts in a list of mutants
        return a dictionary of org_id --> re_id+co1+co2 --> row
        '''
        organism = Organism(self.dbname)
        
        if not organism.isOrg(ref_id):
            logger.warning('Organism %s is not present yet!'%ref_id)
            raise Exception('This Organism (%s) is not present yet!'%ref_id)
        
        muts = set(muts)
        for mut_id in muts:
            if not organism.isOrg(mut_id) or not organism.isMutant(mut_id):
                logger.warning('Wrong organism! (%s)'%mut_id)
                raise Exception('Wrong organism! (%s)'%mut_id)
        
        out = {}
        out[ref_id] = {}
        for row in self.getAllRPairsReacts(ref_id, path_id):
            out[ref_id] = out.get(ref_id, {})
            out[ref_id][row.re_id+row.co1+row.co2] = row
        
        for mut_id in muts:
            # which kind of mutant is this?
            kind = organism.getOrg(mut_id).mkind
            
            if kind == 'insertion':
                ref = out[ref_id]
                mut = {}
                for row in self.getAllRPairsReacts(mut_id, path_id):
                    mut[row.re_id+row.co1+row.co2] = row
                mut = dict(mut.items() + ref.items())
                
                out[mut_id] = {}
                
                keys = set(mut.keys()).difference(set(ref.keys()))
                for k in keys:
                    out[mut_id][k] = mut[k]
            else:
                ref_temp = {}
                for row in self.getAllRPairsReacts(ref_id, path_id):
                    ref_temp[row.re_id+row.co1+row.co2] = row
                
                mut = {}
                for row in self.getAllRPairsReacts(mut_id, path_id):
                    mut[row.re_id+row.co1+row.co2] = row
                    
                for k in mut:
                    if k in ref_temp:
                        del ref_temp[k]
                    
                present = set()
                for k, r in ref_temp.iteritems():
                    present.add(r.re_id+r.co1+r.co2)
                absent = set()
                for k, r in mut.iteritems():
                    absent.add(r.re_id+r.co1+r.co2)
                
                out[mut_id] = {}
                
                for k in absent.difference(present):
                    if k in out[ref_id]:
                        out[mut_id][k] = out[ref_id][k]
                    else:
                        out[mut_id][k] = mut[k]
           
        return out
    
    def getCoreReact(self):
        '''
        Get core genome reactions (and numerosity)
        '''
        nOrg = Organism(self.dbname).howMany()
        genome = Genome(self.dbname)
        
        query = '''
                select distinct re_id, o.group_id
                from ko_react k, mapko m, ortholog o
                where k.ko_id = m.ko_id
                and o.prot_id = m.prot_id;
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query)
                        
        rall = [Row(res, cursor.description) for res in cursor]
        groups = set([x.group_id for x in genome.getCore()])
        
        for r in filter(lambda x: x.group_id in groups, rall):
            setattr(r, 'num', nOrg)
            yield r
            
    def getCoreRPairsReact(self, path_id=None):
        '''
        Get core genome rpairs reactions (and numerosity)
        '''
        nOrg = Organism(self.dbname).howMany()
        genome = Genome(self.dbname)
        
        if not path_id:
            query = '''
                    select distinct k.re_id, co1, co2, o.group_id, re.name
                    from ko_react k, mapko m, protein p, rpair_react rr, rpair rp, reaction re, ortholog o
                    where k.ko_id = m.ko_id
                    and p.prot_id = m.prot_id
                    and p.prot_id = o.prot_id
                    and k.re_id = rr.re_id
                    and rr.rp_id = rp.rp_id
                    and rr.re_id=re.re_id
                    and kind like "%main%";
                    '''
        else:
            query = '''
                    select distinct k.re_id, co1, co2, o.group_id, re.name
                    from ko_react k, mapko m, protein p, rpair_react rr, rpair rp, reaction re, ortholog o, react_path p1
                    where k.ko_id = m.ko_id
                    and p.prot_id = m.prot_id
                    and p.prot_id = o.prot_id
                    and k.re_id = rr.re_id
                    and rr.rp_id = rp.rp_id
                    and rr.re_id=re.re_id
                    and re.re_id=p1.re_id
                    and p1.path_id=?
                    and kind like "%main%";
                    '''
    
        with self.connection as conn:
            if not path_id:
                cursor=conn.execute(query)
            else:
                cursor=conn.execute(query,
                                    [path_id,])
                
        rall = [Row(res, cursor.description) for res in cursor]
        groups = set([x.group_id for x in genome.getCore()])
        
        already = set()
        for r in filter(lambda x: x.group_id in groups, rall):
            setattr(r, 'weight', nOrg)
            if r.re_id+r.co1+r.co2 in already:
                continue
            already.add(r.re_id+r.co1+r.co2)
            yield r
    
    def getDispensableReact(self):
        '''
        Get dispensable genome reactions (and numerosity)
        '''
        genome = Genome(self.dbname)
        
        query = '''
                select distinct re_id, o.group_id
                from ko_react k, mapko m, ortholog o, protein p
                where k.ko_id = m.ko_id
                and o.prot_id = m.prot_id
                and o.prot_id = p.prot_id;
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query)
            
        rall = [Row(res, cursor.description) for res in cursor]
        groups = set([x.group_id for x in genome.getDisp()])
        
        for r in filter(lambda x: x.group_id in groups, rall):
            setattr(r, 'num', genome.getGroupNum(r.group_id))
            yield r
    
    def getDispensableRPairsReact(self, path_id=None):
        '''
        Get dispensable genome rpairs reactions (and numerosity)
        '''
        genome = Genome(self.dbname)
        
        if not path_id:
            query = '''
                select distinct k.re_id, co1, co2, o.group_id, re.name
                from ko_react k, mapko m, protein p, rpair_react rr, rpair rp, reaction re, ortholog o
                where k.ko_id = m.ko_id
                and p.prot_id = m.prot_id
                and p.prot_id = o.prot_id
                and k.re_id = rr.re_id
                and rr.rp_id = rp.rp_id
                and rr.re_id=re.re_id
                and kind like "%main%";
                '''
        else:
            query = '''
                select distinct k.re_id, co1, co2, o.group_id, re.name
                from ko_react k, mapko m, protein p, rpair_react rr, rpair rp, reaction re, ortholog o, react_path p1
                where k.ko_id = m.ko_id
                and p.prot_id = m.prot_id
                and p.prot_id = o.prot_id
                and k.re_id = rr.re_id
                and rr.rp_id = rp.rp_id
                and rr.re_id=re.re_id
                and kind like "%main%"
                and re.re_id = p1.re_id
                and p1.path_id = ?;
                '''
                
        with self.connection as conn:
            if not path_id:
                cursor=conn.execute(query)
            else:
                cursor=conn.execute(query,
                                    [path_id,])
            
        rall = [Row(res, cursor.description) for res in cursor]
        groups = set([x.group_id for x in genome.getDisp()])
        
        # Get max organism numerosity
        rpath = {}
        for r in filter(lambda x: x.group_id in groups, rall):
            rpath[r.re_id+r.co1+r.co2] = rpath.get(r.re_id+r.co1+r.co2, set())
            rpath[r.re_id+r.co1+r.co2].add(genome.getGroupNum(r.group_id))
        
        already = set()
        for r in filter(lambda x: x.group_id in groups, rall):
            if r.re_id+r.co1+r.co2 in already:
                continue
            already.add(r.re_id+r.co1+r.co2)
            setattr(r, 'weight', max([w for w in rpath[r.re_id+r.co1+r.co2]]))
            yield r
    
    def getAccessoryReact(self):
        '''
        Get accessory genome reactions (and numerosity)
        '''
        genome = Genome(self.dbname)
        
        query = '''
                select distinct re_id, o.group_id
                from ko_react k, mapko m, ortholog o, protein p
                where k.ko_id = m.ko_id
                and o.prot_id = m.prot_id
                and o.prot_id = p.prot_id;
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query)
            
        rall = [Row(res, cursor.description) for res in cursor]
        groups = set([x.group_id for x in genome.getAcc()])
        
        for r in filter(lambda x: x.group_id in groups, rall):
            setattr(r, 'num', genome.getGroupNum(r.group_id))
            yield r
    
    def getAccessoryRPairsReact(self, path_id=None):
        '''
        Get accessory genome rpairs reactions (and numerosity)
        '''
        genome = Genome(self.dbname)
        
        if not path_id:
            query = '''
                select distinct k.re_id, co1, co2, o.group_id, re.name
                from ko_react k, mapko m, protein p, rpair_react rr, rpair rp, reaction re, ortholog o
                where k.ko_id = m.ko_id
                and p.prot_id = m.prot_id
                and p.prot_id = o.prot_id
                and k.re_id = rr.re_id
                and rr.rp_id = rp.rp_id
                and rr.re_id=re.re_id
                and kind like "%main%";
                '''
        else:
            query = '''
                select distinct k.re_id, co1, co2, o.group_id, re.name
                from ko_react k, mapko m, protein p, rpair_react rr, rpair rp, reaction re, ortholog o, react_path p1
                where k.ko_id = m.ko_id
                and p.prot_id = m.prot_id
                and p.prot_id = o.prot_id
                and k.re_id = rr.re_id
                and rr.rp_id = rp.rp_id
                and rr.re_id=re.re_id
                and kind like "%main%"
                and re.re_id = p1.re_id
                and p1.path_id = ?;
                '''
                
        with self.connection as conn:
            if not path_id:
                cursor=conn.execute(query)
            else:
                cursor=conn.execute(query,
                                    [path_id,])
            
        rall = [Row(res, cursor.description) for res in cursor]
        groups = set([x.group_id for x in genome.getAcc()])
        
        # Get max organism numerosity
        rpath = {}
        for r in filter(lambda x: x.group_id in groups, rall):
            rpath[r.re_id+r.co1+r.co2] = rpath.get(r.re_id+r.co1+r.co2, set())
            rpath[r.re_id+r.co1+r.co2].add(genome.getGroupNum(r.group_id))
        
        already = set()
        for r in filter(lambda x: x.group_id in groups, rall):
            if r.re_id+r.co1+r.co2 in already:
                continue
            already.add(r.re_id+r.co1+r.co2)
            setattr(r, 'weight', max([w for w in rpath[r.re_id+r.co1+r.co2]]))
            yield r
    
    def getUniqueReact(self):
        '''
        Get unique genome reactions (and numerosity)
        '''
        genome = Genome(self.dbname)
        
        query = '''
                select distinct re_id, o.group_id
                from ko_react k, mapko m, ortholog o
                where k.ko_id = m.ko_id
                and o.prot_id = m.prot_id;
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query)
            
        rall = [Row(res, cursor.description) for res in cursor]
        groups = set([x.group_id for x in genome.getUni()])
        
        for r in filter(lambda x: x.group_id in groups, rall):
            setattr(r, 'num', 1)
            yield r
            
    def getUniqueRPairsReact(self, path_id=None):
        '''
        Get unique genome rpairs reactions (and numerosity)
        '''
        genome = Genome(self.dbname)
        
        if not path_id:
            query = '''
                select distinct k.re_id, co1, co2, o.group_id, re.name
                from ko_react k, mapko m, protein p, rpair_react rr, rpair rp, reaction re, ortholog o
                where k.ko_id = m.ko_id
                and p.prot_id = m.prot_id
                and p.prot_id = o.prot_id
                and k.re_id = rr.re_id
                and rr.rp_id = rp.rp_id
                and rr.re_id=re.re_id
                and kind like "%main%";
                '''
        else:
            query = '''
                select distinct k.re_id, co1, co2, o.group_id, re.name
                from ko_react k, mapko m, protein p, rpair_react rr, rpair rp, reaction re, ortholog o, react_path p1
                where k.ko_id = m.ko_id
                and p.prot_id = m.prot_id
                and p.prot_id = o.prot_id
                and k.re_id = rr.re_id
                and rr.rp_id = rp.rp_id
                and rr.re_id=re.re_id
                and kind like "%main%"
                and re.re_id = p1.re_id
                and p1.path_id = ?;
                '''
                
        with self.connection as conn:
            if not path_id:
                cursor=conn.execute(query)
            else:
                cursor=conn.execute(query,
                                    [path_id,])
            
        rall = [Row(res, cursor.description) for res in cursor]
        groups = set([x.group_id for x in genome.getUni()])
        
        already = set()
        for r in filter(lambda x: x.group_id in groups, rall):
            if r.re_id+r.co1+r.co2 in already:
                continue
            already.add(r.re_id+r.co1+r.co2)
            setattr(r, 'weight', 1)
            yield r
            
    def getOrgReact(self, org_id):
        '''
        Get reactions from a defined organism (and numerosity)
        '''
        query = '''
                select distinct re_id, count(distinct p.prot_id) num
                from ko_react k, mapko m, protein p
                where k.ko_id = m.ko_id
                and p.prot_id = m.prot_id
                and org_id = ?
                group by re_id
                order by num DESC;
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query,[org_id,])
            
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getReactOrg(self, re_id):
        '''
        Get organism(s) from a defined reaction
        '''
        query = '''
                select distinct org_id
                from ko_react k, mapko m, protein p
                where k.ko_id = m.ko_id
                and p.prot_id = m.prot_id
                and re_id = ?
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query,[re_id,])
            
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getReferenceReact(self, mut_id, ref_id):
        '''
        Get reactions from a reference organism (and numerosity)
        The mutated proteins won't be taken into account
        '''
        query = '''
                select distinct re_id, count(distinct p.prot_id) num
                from ko_react k, mapko m, protein p
                where k.ko_id = m.ko_id
                and p.prot_id = m.prot_id
                and org_id = ?
                and p.prot_id not in (select prot_id
                                    from protein
                                    where org_id = ?)
                group by re_id
                order by num DESC;
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query,[ref_id,mut_id,])
            
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getExclusiveReactions(self, orgs=set()):
        '''
        Return the number of reactions ID exclusive to a list of organisms
        return a dictionary of org_id --> set(re_id, ...)
        if the organisms list is empty, all the organisms are queried
        '''
        organism = Organism(self.dbname)
        if len(orgs) == 0:
            orgs = [org.org_id for org in organism.getAll()]
        else:
            orgs = set(orgs)
            for org_id in orgs:
                if not organism.isOrg(org_id):
                    logger.warning('Organism %s is not present yet!'%org_id)
                    raise Exception('This Organism (%s) is not present yet!'%org_id)
        
        react = {}
        
        for org_id in orgs:
            query = '''
                    select distinct re_id
                    from protein p, mapko m, ko_react k
                    where org_id=?
                    and p.prot_id=m.prot_id
                    and m.ko_id=k.ko_id;
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query,[org_id,])
            
                r = set()
                for res in cursor:
                    r.add(res[0])
                react[org_id] = r
        
        out = {}
        
        # Thanks to @balpha for this one
        # http://stackoverflow.com/a/2042394/1237531
        from itertools import combinations
        all_orgs = set(react.keys())
        for i in combinations(react, len(react)-1):
            this_orgs = set(i)
            org_id = all_orgs.difference(this_orgs).pop()
            out[org_id] = react[org_id].difference(*[react[x] for x in this_orgs])            
           
        return out
    
    def getExclusiveReactionsMutants(self, ref_id, muts=set()):
        '''
        Return the reactions ID exclusive to a list of organisms
        return a dictionary of org_id --> set(re_id, ...)
        '''
        organism = Organism(self.dbname)
        
        if not organism.isOrg(ref_id):
            logger.warning('Organism %s is not present yet!'%ref_id)
            raise Exception('This Organism (%s) is not present yet!'%ref_id)
        
        muts = set(muts)
        for mut_id in muts:
            if not organism.isOrg(mut_id) or not organism.isMutant(mut_id):
                logger.warning('Wrong organism! (%s)'%mut_id)
                raise Exception('Wrong organism! (%s)'%mut_id)
        
        out = {}
        out[ref_id] = set()
        
        for mut_id in muts:
            # which kind of mutant is this?
            kind = organism.getOrg(mut_id).mkind
            
            if kind == 'insertion':
                ref = self.getExclusiveReactions( [ref_id] )[ref_id]
                mut = self.getExclusiveReactions( [mut_id] )[mut_id].union(ref)
                out[mut_id] = mut.difference(ref)
            else:
                # prot_id --> [re_id, ...]
                ref_react = {}
                for re in self.getAllReactions(ref_id):
                    ref_react[re.prot_id] = ref_react.get(re.prot_id, set())
                    ref_react[re.prot_id].add(re.re_id)
                
                mut_react = {}
                for re in self.getAllReactions(mut_id):
                    mut_react[re.prot_id] = mut_react.get(re.prot_id, set())
                    mut_react[re.prot_id].add(re.re_id)
                    
                for prot_id in mut_react:
                    if prot_id in ref_react:
                        del ref_react[prot_id]
                    
                present = set()
                for prot_id, re_ids in ref_react.iteritems():
                    for re_id in re_ids:
                        present.add(re_id)
                absent = set()
                for prot_id, re_ids in mut_react.iteritems():
                    for re_id in re_ids:
                        absent.add(re_id)
                
                out[mut_id] = absent.difference(present)       
           
        return out
           
    def getExclusiveReactionsPanGenome(self):
        '''
        Return the reactions ID exclusive to each pangenome category
        core, dispensable, accessory, unique
        note: the dispensable genome includes the accessory and the unique
        '''
        genome = Genome(self.dbname)
        
        query = '''
                select distinct k.re_id, o.group_id
                from mapko m, ortholog o, ko_react k
                where m.prot_id = o.prot_id
                and m.ko_id = k.ko_id;
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query)
        
        core = set([x.group_id for x in genome.getCore()])
        disp = set([x.group_id for x in genome.getDisp()])
        acc = set([x.group_id for x in genome.getAcc()])
        uni = set([x.group_id for x in genome.getUni()])
            
        rall = [Row(res, cursor.description) for res in cursor]

        rcore = set([r.re_id for r in filter(lambda x: x.group_id in core, rall)])
        rdisp = set([r.re_id for r in filter(lambda x: x.group_id in disp, rall)])
        racc = set([r.re_id for r in filter(lambda x: x.group_id in acc, rall)])
        runi = set([r.re_id for r in filter(lambda x: x.group_id in uni, rall)])
        
        return (rcore.difference(rdisp),
                rdisp.difference(rcore),
                racc.difference(rcore, runi),
                runi.difference(rcore, racc))
            
    def howManyMapped(self, org_id=None, pangenome=''):
        '''
        Returns the number of proteins mapped to kegg
        If no org_id is provided, the mapped proteins from all organism
        is returned
        '''
        if org_id:
            query = '''
                select count(distinct  p.prot_id)
                from mapko m, protein p
                where m.prot_id = p.prot_id
                and org_id = ?;
                '''
        elif pangenome in ['core', 'dispensable', 'accessory', 'unique']:
            genome = Genome(self.dbname)
            
            query = '''
                    select distinct o.group_id
                    from mapko m, ortholog o
                    where m.prot_id = o.prot_id;
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query)
            
            rall = [Row(res, cursor.description) for res in cursor]
            if pangenome == 'core':
                groups = set([x.group_id for x in genome.getCore()])
                
            elif pangenome == 'dispensable':
                groups = set([x.group_id for x in genome.getDisp()])
                
            elif pangenome == 'accessory':
                groups = set([x.group_id for x in genome.getAcc()])    

            elif pangenome == 'unique':
                groups = set([x.group_id for x in genome.getUni()])
                
            mapped = set([r.group_id
                         for r in filter(lambda x: x.group_id in groups, rall)])
            return len(mapped)
        
        else:
            query = '''
                    select count(distinct  m.prot_id)
                    from mapko m;
                    '''
            
        with self.connection as conn:
            if org_id:
                cursor=conn.execute(query,[org_id,])
            else:
                cursor=conn.execute(query)
        return int(cursor.fetchall()[0][0])
    
    def howManyKO(self, org_id=None, pangenome=''):
        '''
        Returns the number of KO mapped
        If no org_id is provided, the whole number of KO from all organism
        is returned
        '''
        if org_id:
            query = '''
                select count(distinct ko_id)
                from mapko m, protein p
                where m.prot_id = p.prot_id
                and org_id = ?
                '''
        elif pangenome in ['core', 'dispensable', 'accessory', 'unique']:
            genome = Genome(self.dbname)
            
            query = '''
                    select distinct ko_id, o.group_id
                    from mapko m, ortholog o
                    where m.prot_id = o.prot_id;
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query)
            
            rall = [Row(res, cursor.description) for res in cursor]
            if pangenome == 'core':
                groups = set([x.group_id for x in genome.getCore()])
                
            elif pangenome == 'dispensable':
                groups = set([x.group_id for x in genome.getDisp()])
                
            elif pangenome == 'accessory':
                groups = set([x.group_id for x in genome.getAcc()])    

            elif pangenome == 'unique':
                groups = set([x.group_id for x in genome.getUni()])
                
            kos = set([r.ko_id
                         for r in filter(lambda x: x.group_id in groups, rall)])
            return len(kos)
        else:
            query = '''
                    select count(distinct ko_id)
                    from mapko m
                    '''
            
        with self.connection as conn:
            if org_id:
                cursor=conn.execute(query,[org_id,])
            else:
                cursor=conn.execute(query)
        return int(cursor.fetchall()[0][0])
    
    def howManyReactions(self, org_id=None, pangenome=''):
        '''
        Returns the number of reactions mapped
        If no org_id is provided, the whole number of reactions from all organism
        is returned
        '''
        if org_id:
            query = '''
                select count(k.re_id)
                from mapko m, protein p, ko_react k
                where m.prot_id = p.prot_id
                and org_id = ?
                and m.ko_id = k.ko_id
                '''
        elif pangenome in ['core', 'dispensable', 'accessory', 'unique']:
            genome = Genome(self.dbname)
            
            query = '''
                    select distinct re_id, o.group_id
                    from mapko m, ortholog o, ko_react k
                    where m.prot_id = o.prot_id
                    and m.ko_id = k.ko_id
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query)
            
            rall = [Row(res, cursor.description) for res in cursor]
            if pangenome == 'core':
                groups = set([x.group_id for x in genome.getCore()])
                
            elif pangenome == 'dispensable':
                groups = set([x.group_id for x in genome.getDisp()])
                
            elif pangenome == 'accessory':
                groups = set([x.group_id for x in genome.getAcc()])    

            elif pangenome == 'unique':
                groups = set([x.group_id for x in genome.getUni()])
                
            pairs = set([(r.re_id, r.group_id)
                         for r in filter(lambda x: x.group_id in groups, rall)])
            return len(pairs)    

        else:
            query = '''
                    select count(k.re_id)
                    from mapko m, ko_react k
                    where m.ko_id = k.ko_id
                    '''
            
        with self.connection as conn:
            if org_id:
                cursor=conn.execute(query,[org_id,])
            else:
                cursor=conn.execute(query)
        return int(cursor.fetchall()[0][0])
    
    def howManyUniqueReactions(self, org_id=None, pangenome=''):
        '''
        Returns the number of unique reaction IDs are mapped
        If no org_id is provided, the whole number of reactions from all organism
        is returned
        '''
        if org_id:
            query = '''
                    select count(distinct k.re_id)
                    from mapko m, protein p, ko_react k
                    where m.prot_id = p.prot_id
                    and org_id = ?
                    and m.ko_id = k.ko_id
                    '''
        elif pangenome in ['core', 'dispensable', 'accessory', 'unique']:
            genome = Genome(self.dbname)
            
            query = '''
                    select distinct k.re_id, o.group_id
                    from mapko m, ortholog o, ko_react k
                    where m.prot_id = o.prot_id
                    and m.ko_id = k.ko_id;
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query)
            
            rall = [Row(res, cursor.description) for res in cursor]
            if pangenome == 'core':
                groups = set([x.group_id for x in genome.getCore()])
                
            elif pangenome == 'dispensable':
                groups = set([x.group_id for x in genome.getDisp()])
                
            elif pangenome == 'accessory':
                groups = set([x.group_id for x in genome.getAcc()])    

            elif pangenome == 'unique':
                groups = set([x.group_id for x in genome.getUni()])
                
            reacts = set([r.re_id for r in filter(lambda x: x.group_id in groups, rall)])
            return len(reacts)
            
        else:
            query = '''
                    select count(distinct k.re_id)
                    from mapko m, ko_react k
                    where m.ko_id = k.ko_id
                    '''
            
        with self.connection as conn:
            if org_id:
                cursor=conn.execute(query,[org_id,])
            else:
                cursor=conn.execute(query)
        return int(cursor.fetchall()[0][0])
    
    def howManyPathways(self, org_id=None, pangenome=''):
        '''
        Returns the number of pathways mapped
        If no org_id is provided, the whole number of pathways from all organism
        is returned
        '''
        if org_id:
            query = '''
                select count(distinct  path_id)
                from mapko m, protein p, ko_react k, react_path r
                where m.prot_id = p.prot_id
                and org_id = ?
                and m.ko_id = k.ko_id
                and k.re_id = r.re_id
                '''
        
        elif pangenome in ['core', 'dispensable', 'accessory', 'unique']:
            genome = Genome(self.dbname)
            
            query = '''
                    select distinct  path_id, o.group_id
                    from mapko m, ortholog o, ko_react k, react_path r
                    where m.prot_id = o.prot_id
                    and m.ko_id = k.ko_id
                    and k.re_id = r.re_id;
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query)
            
            rall = [Row(res, cursor.description) for res in cursor]
            if pangenome == 'core':
                groups = set([x.group_id for x in genome.getCore()])
                
            elif pangenome == 'dispensable':
                groups = set([x.group_id for x in genome.getDisp()])
                
            elif pangenome == 'accessory':
                groups = set([x.group_id for x in genome.getAcc()])    

            elif pangenome == 'unique':
                groups = set([x.group_id for x in genome.getUni()])
                
            paths = set([r.path_id for r in filter(lambda x: x.group_id in groups, rall)])
            return len(paths)
        
        else:
            query = '''
                    select count(distinct  path_id)
                    from mapko m, ko_react k, react_path r
                    where m.ko_id = k.ko_id
                    and k.re_id = r.re_id
                    '''
            
        with self.connection as conn:
            if org_id:
                cursor=conn.execute(query,[org_id,])
            else:
                cursor=conn.execute(query)
        return int(cursor.fetchall()[0][0])
    
    def getMappedPathways(self, org_id=None, pangenome=''):
        '''
        Generator to the mapped pathways
        If no org_id is provided, the whole number of pathways from all organism
        is returned
        '''
        if org_id:
            query = '''
                select distinct p.path_id, p.name
                from mapko m, protein p, ko_react k, react_path r, pathway p
                where m.prot_id = p.prot_id
                and org_id = ?
                and m.ko_id = k.ko_id
                and k.re_id = r.re_id
                and r.path_id = p.path_id
                '''
        
        elif pangenome in ['core', 'dispensable', 'accessory', 'unique']:
            genome = Genome(self.dbname)
            
            query = '''
                    select distinct p.path_id, p.name
                    from mapko m, ortholog o, ko_react k, react_path r, pathway p
                    where m.prot_id = o.prot_id
                    and m.ko_id = k.ko_id
                    and k.re_id = r.re_id
                    and r.path_id = p.path_id;
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query)
            
            rall = [Row(res, cursor.description) for res in cursor]
            if pangenome == 'core':
                groups = set([x.group_id for x in genome.getCore()])
                
            elif pangenome == 'dispensable':
                groups = set([x.group_id for x in genome.getDisp()])
                
            elif pangenome == 'accessory':
                groups = set([x.group_id for x in genome.getAcc()])    

            elif pangenome == 'unique':
                groups = set([x.group_id for x in genome.getUni()])
                
            paths = set([r.path_id for r in filter(lambda x: x.group_id in groups, rall)])
            for p in paths:
                yield Row([p], ['path_id'])
        
        else:
            query = '''
                    select distinct p.path_id, p.name
                    from mapko m, ko_react k, react_path r, pathway p
                    where m.ko_id = k.ko_id
                    and k.re_id = r.re_id
                    and r.path_id = p.path_id
                    '''
            
        with self.connection as conn:
            if org_id:
                cursor=conn.execute(query,[org_id,])
            else:
                cursor=conn.execute(query)
        
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getConservedReactions(self, path_id=None):
        '''
        Get the reactions that are present in each organism
        This does not consider the orthologs but just the reaction IDs
        '''
        # How many organisms are present?
        oCheck = Organism(self.dbname)
        nOrgs = oCheck.howMany()
        
        if path_id is None:
            query = '''
                    select distinct re_id, count(distinct org_id) orgs
                    from protein p, mapko m, ko_react k
                    where p.prot_id=m.prot_id
                    and m.ko_id=k.ko_id
                    group by re_id
                    having orgs=?
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query,[nOrgs,])
        else:
            query = '''
                    select distinct k.re_id, count(distinct org_id) orgs
                    from protein p, mapko m, ko_react k, react_path r
                    where p.prot_id=m.prot_id
                    and m.ko_id=k.ko_id
                    and k.re_id=r.re_id
                    and path_id=?
                    group by k.re_id
                    having orgs=?
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query,[path_id, nOrgs,])
            
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getVariableReactions(self, path_id=None):
        '''
        Get the reactions that are differentially present in each organism
        This does not consider the orthologs but just the reaction IDs
        '''
        # How many organisms are present?
        oCheck = Organism(self.dbname)
        nOrgs = oCheck.howMany()
        
        if path_id is None:
            query = '''
                    select distinct re_id, count(distinct org_id) orgs
                    from protein p, mapko m, ko_react k
                    where p.prot_id=m.prot_id
                    and m.ko_id=k.ko_id
                    group by re_id
                    having orgs<?
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query,[nOrgs,])
        else:
            query = '''
                    select distinct k.re_id, count(distinct org_id) orgs
                    from protein p, mapko m, ko_react k, react_path r
                    where p.prot_id=m.prot_id
                    and m.ko_id=k.ko_id
                    and k.re_id=r.re_id
                    and path_id=?
                    group by k.re_id
                    having orgs<?
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query,[path_id, nOrgs,])
            
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getConservedRPairsReact(self, path_id=None):
        '''
        Get the reactions that are present in each organism
        This does not consider the orthologs but just the reaction IDs
        '''
        # How many organisms are present?
        oCheck = Organism(self.dbname)
        nOrgs = oCheck.howMany()
        
        if path_id is None:
            query = '''
                    select k.re_id, co1, co2, re.name, count(distinct org_id) orgs
                    from protein p, mapko m, ko_react k, rpair_react rr, rpair rp, reaction re
                    where p.prot_id=m.prot_id
                    and m.ko_id=k.ko_id
                    and rr.re_id=k.re_id
                    and rr.rp_id=rp.rp_id
                    and rr.re_id=re.re_id
                    group by k.re_id
                    having orgs=?
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query,[nOrgs,])
        else:
            query = '''
                    select k.re_id, co1, co2, re.name, count(distinct org_id) orgs
                    from protein p, mapko m, ko_react k, rpair_react rr, rpair rp, reaction re, react_path r
                    where p.prot_id=m.prot_id
                    and m.ko_id=k.ko_id
                    and rr.re_id=k.re_id
                    and rr.rp_id=rp.rp_id
                    and rr.re_id=re.re_id
                    and k.re_id=r.re_id
                    and path_id=?
                    group by k.re_id
                    having orgs=?
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query,[path_id, nOrgs,])
            
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getVariableRPairsReact(self, path_id=None):
        '''
        Get the reactions that are differentially present in each organism
        This does not consider the orthologs but just the reaction IDs
        '''
        # How many organisms are present?
        oCheck = Organism(self.dbname)
        nOrgs = oCheck.howMany()
        
        if path_id is None:
            query = '''
                    select k.re_id, co1, co2, re.name, count(distinct org_id) orgs
                    from protein p, mapko m, ko_react k, rpair_react rr, rpair rp, reaction re
                    where p.prot_id=m.prot_id
                    and m.ko_id=k.ko_id
                    and rr.re_id=k.re_id
                    and rr.rp_id=rp.rp_id
                    and rr.re_id=re.re_id
                    group by k.re_id
                    having orgs<?
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query,[nOrgs,])
        else:
            query = '''
                    select k.re_id, co1, co2, re.name, count(distinct org_id) orgs
                    from protein p, mapko m, ko_react k, rpair_react rr, rpair rp, reaction re, react_path r
                    where p.prot_id=m.prot_id
                    and m.ko_id=k.ko_id
                    and rr.re_id=k.re_id
                    and rr.rp_id=rp.rp_id
                    and rr.re_id=re.re_id
                    and k.re_id=r.re_id
                    and path_id=?
                    group by k.re_id
                    having orgs<?
                    '''
            
            with self.connection as conn:
                cursor=conn.execute(query,[path_id, nOrgs,])
            
        for res in cursor:
            yield Row(res, cursor.description)
    
class Biolog(DBBase):
    '''
    Class Biolog
    Handles all the data about Biolog entries
    '''
    def __init__(self, dbname='storage'):
        DBBase.__init__(self, dbname)
    
    def resetProject(self):
        '''
        Reset the project Phenomic statuses
        '''
        # Reset the project statuses
        oProj = Project(self.dbname)
        oProj.setPhenome('none')
        
    def clearAllPhenome(self):
        '''
        Truncate all the tables about the phenomic data
        '''
        logger.debug('Clearing phenomic data')
        
        with self.connection as conn:
            conn.execute('delete from biolog_exp;')
            conn.execute('delete from biolog_exp_det;')
            conn.execute('delete from biolog_purged_exp;')
            conn.execute('delete from biolog_purged_exp_det;')
            
        oOrg = Organism(self.dbname)
        oOrg.resetPhenomes()
        
        self.resetProject()
        
    def exportBiolog(self):
        '''
        Generator for biolog data export
        All the relevant data for the biolog plates is extracted
        '''
        with self.connection as conn:
            cursor=conn.execute('select * from biolog;')
        
        yield '# DuctApe generated dump of the Biolog plates'
        yield '#'
        yield '# NOTE: the column "concentration" should be set to 0'
        yield '# set values from 1 to N (i.e. 4) in chemical sensitivity plates'
        yield '#'
        yield '#' + '\t'.join([str(field[0]) 
                             for field in cursor.description])
                             
        for res in cursor:
            yield '\t'.join([str(res[cursor.description.index(field)]) 
                             for field in cursor.description])
                             
    def importBiolog(self, infile):
        '''
        Imports the content of the file object into the biolog table
        Used when first create the database and to add custom plates
        In case of errors there should be a rollback
        '''
        self.boost()
        
        with self.connection as conn:
            for l in open(infile):
                if l.lstrip().startswith('#'):continue
            
                s = l.rstrip('\n').split('\t')
               
                for i in range(len(s)):
                    if s[i] == 'None' or s[i] == '':
                        s[i] = None
                
                values = ''
                for i in range(len(s)):
                    values += '''?, '''
                values = values.rstrip(', ')
                query = '''insert or replace into biolog values (%s);'''%(values)
                
                conn.execute(query, s)
    
    def create(self):
        '''
        Import the standard biolog plates dump
        
        To be used on project init
        '''
        import os
        import ductape.storage.data as data

        fname = os.path.join(os.path.dirname(data.__file__),'biolog.tsv')
        self.importBiolog(fname)

    def isPlate(self, plate_id):
        '''
        Is this plate present?
        '''
        with self.connection as conn:
            cursor=conn.execute('select count(*) from biolog where plate_id=?;',
                                [plate_id,])
        return bool(cursor.fetchall()[0][0])
    
    def isWell(self, well_id):
        '''
        Is this well present?
        '''
        with self.connection as conn:
            cursor=conn.execute('select count(*) from biolog where well_id=?;',
                                [well_id,])
        return bool(cursor.fetchall()[0][0])
    
    def isZeroSubtracted(self, plate_id, well_id, org_id, replica):
        '''
        Is this particular well zero-subtracted?
        '''
        with self.connection as conn:
            cursor=conn.execute('''select zero from biolog_exp 
                                where plate_id=?
                                and well_id=?
                                and org_id=?
                                and replica=?;''',
                                [plate_id,well_id,org_id,replica,])
            
        return bool(cursor.fetchall()[0][0])
    
    def getPlates(self):
        with self.connection as conn:
            cursor=conn.execute('select distinct plate_id from biolog order by plate_id;')
        
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getPlateWells(self, plate_id):
        '''
        Get the (ordered) list of wells from a distinct plate
        '''
        with self.connection as conn:
            cursor=conn.execute('''select distinct well_id
                                   from biolog
                                   where plate_id=?
                                   order by well_id;''',[plate_id,])
        
        for res in cursor:
            yield Row(res, cursor.description).well_id
    
    def getPlate(self, plate_id):
        with self.connection as conn:
            cursor=conn.execute('select * from biolog where plate_id=? order by well_id;',
                                [plate_id,])
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getWells(self):
        with self.connection as conn:
            cursor=conn.execute('''select distinct well_id
                                   from biolog order by well_id;''')
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getWell(self, plate_id, well_id):
        '''
        Get informations on a single well
        '''
        with self.connection as conn:
            cursor=conn.execute('select * from biolog where plate_id=? and well_id=?;',
                                [plate_id,well_id,])
        
        data = cursor.fetchall()
        if len(data) == 0:
            return Row([], cursor.description)
        else:
            return Row(data[0], cursor.description)
        
    def getAllTitles(self):
        '''
        Get the titles for each well
        '''
        with self.connection as conn:
            cursor=conn.execute('''select *
                                from biolog
                                order by plate_id, well_id;''')
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def isMulti(self, plate_id, well_id):
        '''
        Returns True if there are other near wells with the same compound
        '''
        with self.connection as conn:
            cursor=conn.execute('''select concentration
                                from biolog
                                where plate_id=? and well_id=?;''',
                                [plate_id,well_id,])
        
        data = cursor.fetchall()
        if len(data) == 0:
            return False
        else:
            return bool(data[0][0])
        
    def getMulti(self, plate_id, well_id):
        '''
        If the desired well is a multiple one, returns all the related wells
        Otherwise it returns None
        '''
        if not self.isMulti(plate_id, well_id):
            return
        
        mywell = self.getWell(plate_id, well_id)
        
        # Which concentration am i?
        with self.connection as conn:
            cursor=conn.execute('''select * from biolog where plate_id=? 
                                and chemical= ? and
                                (concentration=? or concentration=? 
                                or concentration=? or concentration=?);''',
                                [plate_id, mywell.chemical, 1, 2, 3, 4,])
            
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getCategs(self, active=False):
        '''
        If active is set, only those categories with at least one signal are returned
        '''
        with self.connection as conn:
            if active:
                cursor=conn.execute('''select distinct category
                                    from biolog b, biolog_exp e
                                    where b.plate_id = e.plate_id
                                    order by b.plate_id;''')
            else:
                cursor=conn.execute('''select distinct category
                                    from biolog order by plate_id;''')
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getAll(self):
        with self.connection as conn:
            cursor=conn.execute('''select * from biolog
                                order by plate_id, well_id;''')
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getPlateCategs(self):
        with self.connection as conn:
            cursor=conn.execute('''select distinct plate_id, category
                                from biolog order by plate_id;''')
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getCategByPlate(self, plate_id):
        with self.connection as conn:
            cursor=conn.execute('select category from biolog where plate_id=?;',
                                [plate_id,])
        
        data = cursor.fetchall()
        if len(data) == 0:
            return Row([], cursor.description)
        else:
            return Row(data[0], cursor.description)
        
    def getAllCateg(self, category):
        with self.connection as conn:
            cursor=conn.execute('''select * from biolog where category=?
                                order by plate_id, well_id;''',
                                [category,])
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getByCo(self, co_id):
        with self.connection as conn:
            cursor=conn.execute('''select * from biolog where co_id=?
                                order by plate_id, well_id;''',
                                [co_id,])
        
        for res in cursor:
            yield Row(res, cursor.description)
        
    def getCos(self):
        with self.connection as conn:
            cursor=conn.execute('''select distinct co_id from biolog 
                                where co_id is not null
                                order by co_id;''')
        
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getCosByPlate(self, plate_id):
        with self.connection as conn:
            cursor=conn.execute('''select distinct co_id from biolog 
                                where co_id is not null
                                and plate_id=?
                                order by co_id;''',
                                [plate_id,])
        
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getCosByCateg(self, category):
        with self.connection as conn:
            cursor=conn.execute('''select distinct co_id from biolog 
                                where co_id is not null
                                and category=?
                                order by co_id;''',
                                [category,])
        
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getAllCo(self):
        with self.connection as conn:
            cursor=conn.execute('''select * from biolog 
                                where co_id is not null
                                order by plate_id, well_id;''')
        
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getAllCoByCateg(self, category):
        with self.connection as conn:
            cursor=conn.execute('''select * from biolog 
                                where co_id is not null
                                and category = ?
                                order by plate_id, well_id;''',
                                [category,])
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def addWells(self, explist, clustered=True, replace=False, imported=False):
        '''
        Input: a series of Well objects
        If clustered = True, it is assumed that we have generated the 
        activity parameters and calculated the activity index
        If replace = True, we are merely updating a well
        Checks are performed
        '''
        
        query = '''insert or replace into biolog_exp 
                            (plate_id, well_id, org_id, replica, activity, 
                            zero, min, max, height, plateau, slope, lag,
                            area, v, y0, model, source)
                            values '''
        query1a = '''insert or replace into biolog_exp 
                            (plate_id, well_id, org_id, replica, 
                            zero)
                            values '''
        
        query1 = '''insert or replace into biolog_exp_det
                        (plate_id, well_id, org_id, replica, times, signals)
                        values '''
        
        oCheck = Organism(self.dbname)
        for w in explist:
            if not self.isPlate(w.plate_id):
                logger.warning('Plate %s is not known!'%w.plate_id)
                raise Exception('This plate (%s) is not known!'%w.plate_id)
            if not self.isWell(w.well_id):
                logger.warning('Well %s is not known!'%w.well_id)
                raise Exception('This well (%s) is not known!'%w.well_id)
            if not oCheck.isOrg(w.strain):
                logger.warning('Organism %s is not present yet!'%w.strain)
                raise Exception('This organism (%s) is not present yet!'%w.strain)
            if w.activity is None and clustered and not imported:
                logger.warning('Parameters extraction not yet performed!')
                raise Exception('Parameters extraction not yet performed!')
        
        if not clustered and not replace:
            # Correct the replica
            for w in explist:
                rep = self.howManyReplicasByWell(w.plate_id, w.well_id,
                                                 w.strain)
                w.replica = int(w.replica) + rep
        
        self.boost()
        
        with self.connection as conn:
            if clustered:
                blist = []
                for w in explist:
                    bstr = ('''('%s','%s','%s',%s,'''
                            %(w.plate_id,w.well_id,w.strain,w.replica))
                    for param in [w.activity,int(w.zero),w.min,w.max,
                                  w.height,w.plateau,w.slope,w.lag,
                                  w.area,w.v,w.y0]:
                        if param is None or str(param) == 'nan':
                            bstr += ' null,'
                        else:
                            bstr += ' %s,'%param
                    if w.model is None:
                        bstr += ' null,'
                    else:
                        bstr += ''' '%s','''%w.model
                    if w.source is None:
                        bstr += ' null)'
                    else:
                        bstr += ''' '%s')'''%w.source
                    blist.append(bstr)
            else:
                blist = ['''('%s','%s','%s','%s','%s')'''
                         %(w.plate_id,w.well_id,w.strain,w.replica,int(w.zero))
                         for w in explist]
            
            blist1 = []
            for w in explist:
                blist1 = blist1 + ['''('%s', '%s', '%s', '%s', '%s', '%s')'''
                  %(w.plate_id,w.well_id,w.strain,w.replica,
                   '_'.join([str(x) for x in sorted(w.signals.keys())]),
                   '_'.join([str(w.signals[h]) for h in sorted(w.signals.keys())]))]                
            
            if clustered:
                for bs in get_span(blist, span=1):
                    insert = query + ', '.join(bs)+';'
                    conn.execute(insert)
            else:
                for bs in get_span(blist, span=1):
                    insert = query1a + ', '.join(bs)+';'
                    conn.execute(insert)
                
                for bs in get_span(blist1, span=1):
                    insert = query1 + ', '.join(bs)+';'
                    conn.execute(insert)
            
            conn.execute('''update biolog_exp
                            set model = null where model = '';''')
            conn.execute('''update biolog_exp
                            set source = null where source = '';''')
            conn.commit()
            
    def updateSignals(self, explist):
        '''
        Replaces the signals in the db with the ones provided
        '''
        query = '''update biolog_exp_det
                   set times = ?, signals = ?
                   where plate_id = ?
                   and well_id = ?
                   and org_id = ?
                   and replica = ?'''
        
        self.boost()
        
        with self.connection as conn:
            for w in explist:
                conn.execute(query, 
                             ['_'.join([str(x) for x in sorted(w.signals.keys())]),
                              '_'.join([str(w.signals[h]) for h in sorted(w.signals.keys())]),
                              w.plate_id, w.well_id, w.strain, w.replica])
    
    def delWellsParams(self, wells):
        '''
        Remove all the parameters from the selected wells
        '''
        oCheck = Organism(self.dbname)
        
        query = '''
                update biolog_exp
                set activity=null, 
                    min=null,
                    max=null,
                    height=null,
                    plateau=null,
                    slope=null,
                    lag=null,
                    area=null,
                    v=null,
                    y0=null,
                    model=null,
                    source=null
                where plate_id=? and
                      well_id=? and
                      org_id=? and
                      replica=?;
                '''
        
        for w in wells:
            if not self.isPlate(w.plate_id):
                logger.warning('Plate %s is not known!'%w.plate_id)
                raise Exception('This plate (%s) is not known!'%w.plate_id)
            if not self.isWell(w.well_id):
                logger.warning('Well %s is not known!'%w.well_id)
                raise Exception('This well (%s) is not known!'%w.well_id)
            if not oCheck.isOrg(w.strain):
                logger.warning('Organism %s is not present yet!'%w.strain)
                raise Exception('This organism (%s) is not present yet!'%w.strain)
        
        self.boost()
        
        with self.connection as conn:
            for w in wells:
                conn.execute(query,
                              [w.plate_id,w.well_id,w.strain,w.replica,])
    
    def delWells(self, explist):
        '''
        Input: a series of BiologExp objects
        '''
        with self.connection as conn:
            for w in explist:
                conn.execute('''delete from biolog_exp 
                            where plate_id=? and well_id=? and org_id=?
                            and replica=?;''',
                            [w.plate_id,w.well_id,w.strain,w.replica,])
                
                conn.execute('''delete from biolog_exp_det 
                            where plate_id=? and well_id=? and org_id=?
                            and replica=?;''',
                            [w.plate_id,w.well_id,w.strain,w.replica,])
            conn.commit()
                
    def delOrg(self, org_id):
        '''
        Remove all data about a specific organism
        '''
        with self.connection as conn:
            conn.execute('''delete from biolog_exp 
                        where org_id=?;''',
                        [org_id,])
            
            conn.execute('''delete from biolog_exp_det 
                        where org_id=?;''',
                        [org_id,])
            
            conn.execute('''delete from biolog_purged_exp 
                        where org_id=?;''',
                        [org_id,])
            
            conn.execute('''delete from biolog_purged_exp_det 
                        where org_id=?;''',
                        [org_id,])
        
        org = Organism(self.dbname)
        org.setPhenomeStatus(org_id, 'none')    
    
    def getActivityDistribution(self, plate_id=None):
        '''
        Get the total activity distribution (not replica-wise)
        Returns a dictionary of activity --> #wells
        If plate_id is set, only a particular plate is returned
        '''
        with self.connection as conn:
            if plate_id is None:
                cursor=conn.execute('''select activity, count(*) howmany
                                       from biolog_exp
                                       group by activity
                                       order by activity;''')
            else:
                cursor=conn.execute('''select activity, count(*) howmany
                                   from biolog_exp
                                   where plate_id=?
                                   group by activity
                                   order by activity;''',[plate_id,])
        
        act = {}
        for res in cursor:
            a = Row(res, cursor.description)
            act[a.activity] = a.howmany
            
        return act
    
    def getActivityDistributionByOrg(self, org_id):
        '''
        Get the total activity distribution (not replica-wise)
        Returns a dictionary of activity --> #wells
        '''
        with self.connection as conn:
            cursor=conn.execute('''select activity, count(*) howmany
                                   from biolog_exp
                                   where org_id = ?
                                   group by activity
                                   order by activity;''',
                                   [org_id,])
        
        act = {}
        for res in cursor:
            a = Row(res, cursor.description)
            act[a.activity] = a.howmany
            
        return act
    
    def getActivityDistributionByZero(self, nonzero=False):
        '''
        Get the total activity distribution (not replica-wise)
        if nonzero is set to True, the nonzero subtractable activities are returned
        Returns a dictionary of activity --> #wells
        '''
        with self.connection as conn:
            if not nonzero:
                cursor=conn.execute('''select activity, count(*) howmany
                                       from biolog_exp
                                       where (plate_id='PM01'
                                            or plate_id='PM02A'
                                            or plate_id='PM03B'
                                            or plate_id='PM04A'
                                            or plate_id='PM05'
                                            or plate_id='PM06'
                                            or plate_id='PM07'
                                            or plate_id='PM08'
                                            or plate_id='PM03B'
                                            or plate_id='PM04A')
                                       group by activity
                                       order by activity;''')
            else:
                cursor=conn.execute('''select activity, count(*) howmany
                                       from biolog_exp
                                       where (plate_id!='PM01'
                                            and plate_id!='PM02A'
                                            and plate_id!='PM03B'
                                            and plate_id!='PM04A'
                                            and plate_id!='PM05'
                                            and plate_id!='PM06'
                                            and plate_id!='PM07'
                                            and plate_id!='PM08'
                                            and plate_id!='PM03B'
                                            and plate_id!='PM04A')
                                       group by activity
                                       order by activity;''')
        
        act = {}
        for res in cursor:
            a = Row(res, cursor.description)
            act[a.activity] = a.howmany
            
        return act
    
    def getActivityDistributionByZeroAndOrg(self, org_id, nonzero=False):
        '''
        Get the total activity distribution (not replica-wise)
        if nonzero is set to True, the nonzero subtractable activities are returned
        Returns a dictionary of activity --> #wells
        '''
        with self.connection as conn:
            if not nonzero:
                cursor=conn.execute('''select activity, count(*) howmany
                                       from biolog_exp
                                       where (plate_id='PM01'
                                            or plate_id='PM02A'
                                            or plate_id='PM03B'
                                            or plate_id='PM04A'
                                            or plate_id='PM05'
                                            or plate_id='PM06'
                                            or plate_id='PM07'
                                            or plate_id='PM08'
                                            or plate_id='PM03B'
                                            or plate_id='PM04A')
                                       and org_id = ?
                                       group by activity
                                       order by activity;''',
                                       [org_id,])
            else:
                cursor=conn.execute('''select activity, count(*) howmany
                                       from biolog_exp
                                       where (plate_id!='PM01'
                                            and plate_id!='PM02A'
                                            and plate_id!='PM03B'
                                            and plate_id!='PM04A'
                                            and plate_id!='PM05'
                                            and plate_id!='PM06'
                                            and plate_id!='PM07'
                                            and plate_id!='PM08'
                                            and plate_id!='PM03B'
                                            and plate_id!='PM04A')
                                       and org_id = ?
                                       group by activity
                                       order by activity;''',
                                       [org_id,])
        
        act = {}
        for res in cursor:
            a = Row(res, cursor.description)
            act[a.activity] = a.howmany
            
        return act
    
    def getActivityDistributionByCateg(self, categ):
        '''
        Get the total activity distribution (not replica-wise)
        Returns a dictionary of activity --> #wells
        '''
        with self.connection as conn:
            cursor=conn.execute('''select activity, count(*) howmany
                                   from biolog_exp b1, biolog b
                                   where b.plate_id = b1.plate_id
                                   and b.well_id = b1.well_id
                                   and category = ?
                                   group by activity
                                   order by activity;''',
                                   [categ,])
        
        act = {}
        for res in cursor:
            a = Row(res, cursor.description)
            act[a.activity] = a.howmany
            
        return act
    
    def getActivityDistributionByCategAndOrg(self, categ, org_id):
        '''
        Get the total activity distribution (not replica-wise)
        Returns a dictionary of activity --> #wells
        '''
        with self.connection as conn:
            cursor=conn.execute('''select activity, count(*) howmany
                                   from biolog_exp b1, biolog b
                                   where b.plate_id = b1.plate_id
                                   and b.well_id = b1.well_id
                                   and category = ?
                                   and org_id = ?
                                   group by activity
                                   order by activity;''',
                                   [categ,org_id,])
        
        act = {}
        for res in cursor:
            a = Row(res, cursor.description)
            act[a.activity] = a.howmany
            
        return act
        
    def getMaxActivity(self):
        '''
        Returns the maximum activity present
        '''
        with self.connection as conn:
            cursor=conn.execute('''select max(activity)
                                    from biolog_exp;''')
        
        return cursor.fetchall()[0][0]
    
    def getAvgActivity(self, plate_id, well_id, org_id):
        '''
        Get the average activity for a particular experiment
        '''
        with self.connection as conn:
            cursor=conn.execute('''select avg(activity) from biolog_exp  
                                   where org_id=?
                                   and plate_id=?
                                   and well_id=?;''',
                                  [org_id,plate_id, well_id,])
        
        try:
            return float(cursor.fetchall()[0][0])      
        except:
            return None
    
    def getAvgActivityEachOrg(self, plate_id, well_id):
        '''
        Get the average activity for a particular experiment
        Returns one record for each organism
        '''
        with self.connection as conn:
            cursor=conn.execute('''select org_id, avg(activity) avgact
                                   from biolog_exp  
                                   where plate_id=?
                                   and well_id=?
                                   group by org_id;''',
                                  [plate_id, well_id,])
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getReplicas(self, plate_id, well_id, org_id):
        with self.connection as conn:
            cursor=conn.execute('''select * from biolog_exp  
                        where plate_id=? and well_id=? and org_id=?
                        order by replica;''',
                        [plate_id,well_id,org_id,])
        
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getActiveByPlate(self, plate_id, activity):
        '''
        Get those wells at least active activity
        '''
        with self.connection as conn:
            cursor=conn.execute('''select * from biolog_exp  
                        where plate_id=?
                        and activity>=?;''',
                        [plate_id,activity,])
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getAllActive(self, activity):
        '''
        Get those wells at least active activity
        '''
        with self.connection as conn:
            cursor=conn.execute('''select * from biolog_exp  
                        where activity>=?;''',
                        [activity,])
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def howManyActive(self,activity):
        '''
        How many wells are active at least activity?
        '''
        with self.connection as conn:
            cursor=conn.execute('''select count(*) from biolog_exp
                                    where activity>=?;''',[activity,])
        return int(cursor.fetchall()[0][0])
    
    def howManyActiveByPlate(self, plate_id, activity):
        '''
        How many wells are active at least activity?
        '''
        with self.connection as conn:
            cursor=conn.execute('''select count(*) from biolog_exp
                                    where activity>=?
                                    and plate_id=?;''',[activity,plate_id,])
        return int(cursor.fetchall()[0][0])
    
    def howManyActiveByOrg(self, org_id, activity):
        '''
        How many wells are active at least activity?
        '''
        with self.connection as conn:
            cursor=conn.execute('''select count(*) from biolog_exp
                                    where activity>=?
                                    and org_id=?;''',[activity,org_id,])
        return int(cursor.fetchall()[0][0])
    
    def howManyReplicas(self):
        '''
        How many replicas do we have?
        '''
        with self.connection as conn:
                    cursor=conn.execute('''select count(distinct replica)
                                           from biolog_exp;''')
        return int(cursor.fetchall()[0][0])        
    
    def howManyReplicasByOrg(self, org_id):
        '''
        How many replicas do we have for my organism?
        '''
        with self.connection as conn:
            cursor=conn.execute('''select count(distinct replica)
                                   from biolog_exp
                                   where org_id=?;''',[org_id,])

        return int(cursor.fetchall()[0][0])
    
    def howManyReplicasByWell(self, plate_id, well_id, org_id):
        '''
        How many replicas do we have for this single well?
        '''
        with self.connection as conn:
                    cursor=conn.execute('''select count(distinct replica)
                                           from biolog_exp
                                           where plate_id=?
                                           and well_id=?
                                           and org_id=?;''',[plate_id,
                                                             well_id,
                                                             org_id,])
        return int(cursor.fetchall()[0][0])    
    
    def getOneWell(self, plate_id, well_id, org_id, replica):
        '''
        Get a precise well
        '''
        with self.connection as conn:
            cursor=conn.execute('''select *
                                   from biolog_exp
                                   where org_id=?
                                   and plate_id=?
                                   and well_id=?
                                   and replica=?;''',[org_id,
                                                      plate_id,
                                                      well_id,
                                                      replica,])
        
        data = cursor.fetchall()
        if len(data) == 0:
            return Row([], cursor.description)
        else:
            return Row(data[0], cursor.description)
        
    def getOrgsByWell(self, plate_id, well_id, replica=None):
        '''
        Get the organisms for a precise well
        If replica is provided, only the org_id for that particular replica
        are provided
        '''
        with self.connection as conn:
            if not replica:
                cursor=conn.execute('''select distinct org_id
                                       from biolog_exp
                                       where plate_id=?
                                       and well_id=?;''',[plate_id,
                                                      well_id,])
            else:
                cursor=conn.execute('''select distinct org_id
                                       from biolog_exp
                                       where plate_id=?
                                       and well_id=?
                                       and replica=?;''',[plate_id,
                                                   well_id,
                                                   replica,])                
            
        return [x[0] for x in cursor.fetchall()]
    
    def getAllWells(self):
        '''
        Get all the wells from the storage
        '''
        with self.connection as conn:
            cursor=conn.execute('''select *
                                   from biolog_exp
                                   order by plate_id, well_id, org_id, replica;''')
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getDistinctWells(self, replica=False):
        '''
        Get the distinct wells identifiers (if replica=True, replica aware)
        '''
        with self.connection as conn:
            if replica:
                cursor=conn.execute('''select distinct plate_id, well_id, replica
                                from biolog_exp
                                order by plate_id, well_id, replica;''')
            else:
                cursor=conn.execute('''select distinct plate_id, well_id
                                        from biolog_exp
                                        order by plate_id, well_id;''')                
        
        for res in cursor:
            yield Row(res, cursor.description)        
    
    def getOrgDistinctWells(self, org_id):
        '''
        Get the distinct wells identifiers for a particular organism
        '''
        with self.connection as conn:
            cursor=conn.execute('''select distinct plate_id, well_id
                                   from biolog_exp
                                   where org_id=?
                                   order by plate_id, well_id;''',[org_id,])
        
        for res in cursor:
            yield Row(res, cursor.description)        
    
    def getOrgWells(self, org_id):
        '''
        Get all the wells for a particular organism
        '''
        with self.connection as conn:
            cursor=conn.execute('''select *
                                   from biolog_exp b1
                                   where org_id=?
                                   order by plate_id, well_id,
                                   replica;''',[org_id,])
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def maxSignal(self):
        '''
        Get the maximum signal value
        '''
        with self.connection as conn:
            cursor=conn.execute('select max(max) from biolog_exp;')
        return int(cursor.fetchall()[0][0])
    
    def getRandomWells(self, activity=None, howmany=10, zero=False):
        '''
        Get a random number of wells
        
        May be filtered by activity and zero subtraction state.
        '''
        with self.connection as conn:
            if zero:
                izero = 1
            else:
                izero = 0
            
            if activity is None:
                cursor=conn.execute('''select * from biolog_exp
                                        where zero = ?
                                        order by random()
                                        limit 0,?''', [izero, howmany,])
            else:
                cursor=conn.execute('''select * from biolog_exp
                                        where activity = ?
                                        and zero = ?
                                        order by random()
                                        limit 0,?''',
                                        [activity, izero, howmany,])
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getAllActivity(self):
        '''
        Get all the activities from the storage
        '''
        with self.connection as conn:
            cursor=conn.execute('''select b.plate_id, b.well_id, b.org_id,
                                          b.replica, b.activity
                                   from biolog_exp b;''')
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getAllSignals(self):
        '''
        Get all the signals from the storage
        '''
        with self.connection as conn:
            cursor=conn.execute('''select b.plate_id, b.well_id, b.org_id,
                                          b.replica, b1.times, b1.signals,
                                          b.activity, b.min, b.max, b.height,
                                          b.plateau, b.slope, b.lag, b.area,
                                          b.v, b.y0,
                                          b.model, b.source
                                   from biolog_exp_det b1, biolog_exp b
                                   where b.plate_id=b1.plate_id
                                   and b.well_id=b1.well_id
                                   and b.org_id=b1.org_id
                                   and b.replica=b1.replica;''')
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getAllSignalsNoParams(self):
        '''
        Get all the signals for which we have no parameters from the storage
        '''
        with self.connection as conn:
            cursor=conn.execute('''select b.plate_id, b.well_id, b.org_id,
                                          b.replica, b1.times, b1.signals,
                                          b.activity
                                   from biolog_exp_det b1, biolog_exp b
                                   where b.plate_id=b1.plate_id
                                   and b.well_id=b1.well_id
                                   and b.org_id=b1.org_id
                                   and b.replica=b1.replica
                                   and (activity=null and min=null
                                       and max=null and height=null
                                       and plateau=null and slope=null
                                       and lag=null and area=null
                                       and v=null and y0=null
                                       and model=null and source=null);''')
        
        for res in cursor:
            yield Row(res, cursor.description)
      
    def getParamsSources(self):
        '''
        Generator to the distinct sources for the parameters calculation
        '''
        query = '''
                select distinct source
                from biolog_exp;
                '''
        with self.connection as conn:
            cursor=conn.execute(query)
        
        for res in cursor:
            yield Row(res, cursor.description)
    
    def getZeroSubtractablePlates(self):
        '''
        Generator to the plate IDs that can be zero subtracted
        '''
        query='''select distinct plate_id
                from biolog
                where zero_well_id is not null;
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query)
            
        for res in cursor:
            yield Row(res, cursor.description)
            
    def getControlWells(self):
        '''
        Generator to tuples of plate_id, zero_well
        '''
        query = '''
                select distinct plate_id, well_id
                from biolog
                where zero_well_id is null
                and plate_id in (select distinct plate_id
                                   from biolog
                                   where zero_well_id is not null);
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query)
            
        for res in cursor:
            r = Row(res, cursor.description)
            yield r.plate_id, r.well_id
            
    def getControlPairs(self):
        '''
        Generator to a tuple with plate_id, well_id, zero_well_id
        '''
        query = '''
                select distinct plate_id, well_id, zero_well_id
                from biolog
                where zero_well_id is not null
                and plate_id in (select distinct plate_id
                                   from biolog
                                   where zero_well_id is not null);
                '''
        
        with self.connection as conn:
            cursor=conn.execute(query)
            
        for res in cursor:
            r = Row(res, cursor.description)
            yield r.plate_id, r.well_id, r.zero_well_id
    
    def getZeroSubtractableSignals(self):
        '''
        Get all the signals that can be zero-subtracted
        '''
        with self.connection as conn:
            cursor=conn.execute('''select * from biolog_exp where zero = 0;''')
        
        notYet = [Row(res, cursor.description) for res in cursor]
        for well in notYet:
            with self.connection as conn:
                cursor=conn.execute('''select * from biolog_exp_det
                                    where plate_id=? and well_id=? and org_id=?
                                    and replica=?;''',
                                    [well.plate_id,well.well_id,
                                     well.org_id,well.replica,])
            for res in cursor:
                yield Row(res, cursor.description)
                
    def atLeastOneZeroSubtracted(self):
        '''
        Is there at least one well that is zero subtracted?
        '''
        with self.connection as conn:
            cursor=conn.execute('select count(*) from biolog_exp where zero=1;')
        return bool(cursor.fetchall()[0][0])
    
    def atLeastOneParameter(self):
        '''
        Is there at least one well with the parameters already calculated?
        '''
        with self.connection as conn:
            cursor=conn.execute('''select count(*)
                                    from biolog_exp
                                    where activity is not null;''')
        return bool(cursor.fetchall()[0][0])
    
    def atLeastOneNoParameter(self):
        '''
        Is there at least one well with no parameters already calculated?
        '''
        with self.connection as conn:
            cursor=conn.execute('''select count(*)
                                    from biolog_exp
                                    where activity is null;''')
        return bool(cursor.fetchall()[0][0])
    
    def isEmpty(self):
        '''
        Checks if there is at least some phenomic data
        '''
        with self.connection as conn:
            cursor=conn.execute('''select count(*)
                                    from biolog_exp;''')
        return not bool(cursor.fetchall()[0][0])
    
    def getCompounds2Analyse(self):
        '''
        Get all the biolog compounds ID yet to be anlysed
        '''
        with self.connection as conn:
            cursor=conn.execute('''select distinct co_id from biolog
                                 where co_id is not null
                                 and co_id not in (select co_id
                                                 from compound);''')
            
        for res in cursor:
            yield Row(res, cursor.description)
    
    def atLeastOnePurged(self):
        '''
        Is there at least one well discarded?
        '''
        with self.connection as conn:
            cursor=conn.execute('''select count(*)
                                    from biolog_purged_exp;''')
        return bool(cursor.fetchall()[0][0])
    
    def howManyPurged(self):
        '''
        How many discarded wells?
        '''
        with self.connection as conn:
            cursor=conn.execute('''select count(*)
                                    from biolog_purged_exp;''')
        return int(cursor.fetchall()[0][0])
       
    def moveDiscardedWells(self, wells):
        '''
        Get a list of biolog_ids and move them to the
        "purged wells" zone
        '''
        self.boost()
        
        with self.connection as conn:
            for w in wells:
                cursor = conn.execute('''select * from biolog_exp_det
                                where plate_id=? and well_id=? and org_id=?
                                    and replica=?;''',
                                    [w[0],w[1],w[2],w[3],])
                
                well = Row(cursor.fetchall()[0], cursor.description)
                conn.execute('''insert into biolog_purged_exp_det
                        values (?,?,?,?,?,?);''',[w[0],w[1],w[2],w[3],
                                                  well.times, well.signals])
                conn.execute('''delete from biolog_exp_det
                        where plate_id=? and well_id=? and org_id=?
                        and replica=?;''',
                            [w[0],w[1],w[2],w[3],])
                    
                p, w, o, r = w
                cursor = conn.execute('''select * from biolog_exp
                                where plate_id = ?
                                and well_id = ?
                                and org_id = ?
                                and replica = ?;''',[p,w,o,r,])
                well = Row(cursor.fetchall()[0], cursor.description)
                conn.execute('''insert into biolog_purged_exp
                        values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);''',
                        [p,w,o,r,well.activity, well.zero,well.min,
                         well.max,well.height,well.plateau,
                         well.slope,well.lag,well.area,
                         well.v,well.y0,well.model,well.source])
                conn.execute('''delete from biolog_exp
                                where plate_id = ?
                                and well_id = ?
                                and org_id = ?
                                and replica = ?;''',[p,w,o,r,])
        
    def restoreDiscardedWells(self, plates=[], replica=None):
        '''
        Restore all the discarded wells
        '''
        import copy
        
        self.boost()
        
        restored = 0
        
        with self.connection as conn:
            if replica is None:
                cursor = conn.execute('''select * from biolog_purged_exp_det;''')
            else:
                cursor = conn.execute('''select * from biolog_purged_exp_det
                                        where replica=?;''',[replica,])
                
            exp_det = copy.deepcopy(cursor.description)
            
            for res in cursor:
                well = Row(res, exp_det)
                
                if len(plates) > 0:
                    if well.plate_id not in plates:
                        logger.debug('Skipping restoring of %s'%well.plate_id)
                        continue
                    
                restored += 1
                
                conn.execute('''insert into biolog_exp_det
                        values (?,?,?,?,?,?);''',[well.plate_id, well.well_id,
                                                  well.org_id, well.replica,
                                                  well.times, well.signals])
                conn.execute('''delete from biolog_purged_exp_det
                        where plate_id=? and well_id=? and org_id=?
                        and replica=?;''',
                            [well.plate_id, well.well_id,
                              well.org_id, well.replica,])
                    
                cursor = conn.execute('''select * from biolog_purged_exp
                                where plate_id = ?
                                and well_id = ?
                                and org_id = ?
                                and replica = ?;''',[well.plate_id, well.well_id,
                                                  well.org_id, well.replica,])
                well = Row(cursor.fetchall()[0], cursor.description)
                conn.execute('''insert into biolog_exp
                        values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);''',
                        [well.plate_id, well.well_id,
                         well.org_id, well.replica,
                         well.activity, well.zero,well.min,
                         well.max,well.height,well.plateau,
                         well.slope,well.lag,well.area,
                         well.v,well.y0,well.model,well.source])
                conn.execute('''delete from biolog_purged_exp
                                where plate_id = ?
                                and well_id = ?
                                and org_id = ?
                                and replica = ?;''',[well.plate_id, well.well_id,
                                                  well.org_id, well.replica,])
        return restored
