#!/usr/bin/env python
"""
Blast

Genome library

Classes to handle Blast analysis against a local database
"""

__author__ = "Marco Galardini"

################################################################################
# Imports

import logging
import sys
import subprocess

################################################################################
# Log setup

logger = logging.getLogger('Blast')

################################################################################
# Classes

# Useful class for parsing
class BlastHit:
    def __init__(self,query,align,hsp):
        '''
        Query, Alignment and Hsp are all Biopython objects derived from
        Blast results parsing
        '''
        self.query = query.query
        self.query_id = query.query.split(' ')[0]
        self.query_len = int(query.query_length)
        self.hit = align.hit_id
        self.hit_desc = align.hit_def
        self.hit_len = int(align.length)
        self.identity = float(hsp.identities) / float(hsp.align_length)
        self.align_len = int(hsp.align_length)
        self.mismatches = int(hsp.align_length - hsp.identities - hsp.gaps)
        self.gaps = int(hsp.gaps)
        self.query_start = int(hsp.query_start)
        self.query_end = int(hsp.query_end)
        self.subjct_start = int(hsp.sbjct_start)
        self.subjct_end = int(hsp.sbjct_end)
        self.evalue = float(hsp.expect)
        self.bits = float(hsp.bits)
        
    def getHomologyIndex(self):
        '''
        Get an Index useful for stating the quality of the homology measure
        '''
        import math
        HI=( (math.pow(self.identity,2)*(float(self.hit_len)) /
            (float(self.query_len))*(float(self.align_len)/float(self.query_len)))
            )
        return HI
    
    def getHitCoverage(self):
        '''
        Get the hit coverage
        '''
        return float(float(self.align_len)/float(self.hit_len))
    
    def getQueryCoverage(self):
        '''
        Get the query coverage
        '''
        return float(float(self.align_len)/float(self.query_len))
    
    def getKO(self):
        '''
        Assuming that this hit derives from a KEGG DB
        Returns the KO ID
        '''
        import re
        a=re.search("K[0-9]{1,}",
                    self.hit_desc)
        if a is not None:
            return a.group()
        else:
            return None
        
class Blaster(object):
    def __init__(self):
        self._hits = None
        self._out = ''
        
    def createDB(self,seqFile,dbType,outFile='BlastDB',parseIDs=True,
                        title='Generic Blast DB'):
        '''Generation of a Blast DB'''
        cmd = ('makeblastdb -in %s -dbtype %s -out %s -title "%s"')
        cmd = cmd%(seqFile,dbType,outFile,title)
        if parseIDs:
            cmd = cmd+' -parse_seqids'
        logger.info('Create Blast DB cmd: %s'%cmd)
        proc = subprocess.Popen(cmd,shell=(sys.platform!="win32"),
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
        out = proc.communicate()
        return_code = proc.returncode
        if return_code != 0:
            logger.warning('Blast DB creation failed with error %d'
                            %return_code)
        else:
            logger.debug('Blast DB creation successful')
        return bool(not return_code)
    
    def retrieveFromDB(self, db, accession, out='out.fsa', isFile=False):
        '''Retrieve the desired sequence(s) from a Blast DB'''
        if not isFile:
            cmd=('blastdbcmd -db %s -entry "%s" > %s'
                 %(db,accession,out))
        else:
            cmd=('blastdbcmd -db %s -entry_batch "%s" > %s'
                 %(db,accession,out))
        logger.info('BlastDBcmd cmd: %s'%cmd)
        proc = subprocess.Popen(cmd,shell=(sys.platform!="win32"),
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
        out = proc.communicate()
        return_code = proc.returncode
        if return_code != 0:
            logger.warning('BlastDBcmd failed with error %d'
                            %return_code)
        else:
            logger.debug('BlastDBcmd successful')
        return bool(not return_code)
    
    def runBlast(self, queryFile, db, outFile, evalue = 10,
                    task = '', ncpus = 1, additional = ''):
        '''Run Blast with the desired parameters'''
        # Create the command line
        from Bio.Blast.Applications import NcbiblastpCommandline
        self._out = outFile
        cmd = NcbiblastpCommandline(query=queryFile, db=db,
                evalue=float(evalue),
                outfmt='5',out=outFile,
                num_threads=ncpus)
        if task != '':
            cmd.set_parameter('task', task)
        if additional !='':
            cmd = str(cmd)+' '+additional
        cmd=str(cmd)
        logger.info('Run Blast cmd: %s'%cmd)
        # Run Blast and check the return code
        proc = subprocess.Popen(cmd,shell=(sys.platform!="win32"),
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
        out = proc.communicate()
        return_code = proc.returncode
        if return_code != 0:
            logger.warning('Run Blast failed with error %d'
                            %return_code)
        else:
            logger.debug('Run Blast successful')
        return bool(not return_code)
    
    def parseBlast(self, fileOut):
        '''Parse the xml blast output -- default file is self._out'''
        from Bio.Blast import NCBIXML
        self._out = fileOut
        handle = open(fileOut)
        self._hits = NCBIXML.parse(handle)
        
    def getHits(self,expect=10.0):
        '''Returns a Generator query -> BlastObj'''
        if self._hits == None:
            self.parseBlast(self._out)
        for BlastQuery in self._hits:
            hits = []
            for alignment in BlastQuery.alignments:
                for hsp in alignment.hsps:
                    if float(hsp.expect) > expect:continue
                    # Save the hit details
                    h=BlastHit(BlastQuery,alignment,hsp)
                    hits.append(h)
            yield hits
