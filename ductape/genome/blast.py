#!/usr/bin/env python
"""
Blast

Genome library

Classes to handle Blast analysis against a local database
"""
import logging
import os
import subprocess
import sys
try:
    from StringIO import StringIO # Python 2
except ImportError:
    from io import StringIO # Python 3

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('ductape.blast')

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
    def __init__(self, useDisk=False):
        self._hits = None
        self._out = ''
        
        # No-disk
        self._useDisk = bool(useDisk)
        self.retrieved = ''
        self.query = ''
        self.out = ''
        
    def createDB(self,seqFile,dbType,outFile='BlastDB',parseIDs=True,
                        title='Generic Blast DB'):
        '''Generation of a Blast DB'''
        cmd = ('makeblastdb -in %s -dbtype %s -out %s -title "%s"')
        cmd = cmd%(seqFile,dbType,outFile,title)
        if parseIDs:
            cmd = cmd+' -parse_seqids'
        logger.debug('Create Blast DB cmd: %s'%cmd)
        proc = subprocess.Popen(cmd,shell=(sys.platform!="win32"),
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
        out = proc.communicate()
        return_code = proc.returncode
        if return_code != 0:
            logger.warning('Blast DB creation failed with error %d'
                            %return_code)
            logger.warning('%s'%str(out[1]))

        return bool(not return_code)
    
    def retrieveFromDB(self, db, accession, out='out.fsa', isFile=False):
        '''Retrieve the desired sequence(s) from a Blast DB'''
        if not isFile:
            cmd=('blastdbcmd -db %s -entry "%s"'
                 %(db,accession))
        else:
            cmd=('blastdbcmd -db %s -entry_batch "%s"'
                 %(db,accession))
        
        if self._useDisk:
            cmd += ' > %s'%out
        
        logger.debug('BlastDBcmd cmd: %s'%cmd)
        proc = subprocess.Popen(cmd,shell=(sys.platform!="win32"),
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
        out = proc.communicate()
        
        if not self._useDisk:
            self.retrieved = out[0]
        
        return_code = proc.returncode
        if return_code != 0:
            logger.warning('BlastDBcmd failed with error %d'
                            %return_code)
            logger.warning('%s'%str(out[1]))

        return bool(not return_code)
    
    def runBlast(self, queryFile, db, outFile='', evalue = 10,
                    task = '', ncpus = 1, additional = '', outfmt='5'):
        '''Run Blast with the desired parameters'''
        # Create the command line
        from Bio.Blast.Applications import NcbiblastpCommandline
        self._out = outFile
        cmd = NcbiblastpCommandline(db=db,
                evalue=float(evalue),
                outfmt=outfmt,
                num_threads=ncpus)
        if self._useDisk:
            cmd.set_parameter('query', queryFile)
            if outFile != '':
                cmd.set_parameter('out', outFile)
        if task != '':
            cmd.set_parameter('task', task)
        if additional !='':
            cmd = str(cmd)+' '+additional
        cmd=str(cmd)
        logger.debug('Run Blast cmd: %s'%cmd)
        # Run Blast and check the return code
        proc = subprocess.Popen(cmd,shell=(sys.platform!="win32"),
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
        
        if not self._useDisk:
            proc.stdin.write(self.query)
                    
        out = proc.communicate()
        
        if not self._useDisk:
            self.out = out[0]
        
        return_code = proc.returncode
        if return_code != 0:
            logger.warning('Run Blast failed with error %d'
                            %return_code)
            logger.warning('%s'%str(out[1]))

        return bool(not return_code)
    
    def parseBlast(self, fileOut):
        '''Parse the xml blast output -- default file is self._out'''
        from Bio.Blast import NCBIXML

        if self._useDisk:
            self._out = fileOut
            handle = open(fileOut)
        else:
            handle = StringIO(self.out)
            
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
            
class RunBBH(object):
    def __init__(self, query, queryid,
                 source, target, targetorg,
                 evalue, matrix, short = False, uniqueid = 1,
                 kegg = False, ko_entry = None, ko_id = None, useDisk=True):
        self.query = query
        self.queryid = queryid
        self.source = source
        self.target = target
        self.targetorg = targetorg
        self.evalue = evalue
        self.matrix = matrix
        self.short = short
        self.uniqueid = uniqueid
        self.kegg = kegg
        self.ko_entry = ko_entry
        self.ko_id = ko_id
        self.useDisk = bool(useDisk)
        
        self.out = self.query + '_' + str(self.uniqueid) +'.xml'
        self.blaster = Blaster(useDisk=self.useDisk)
        self.additional = (' -soft_masking true -dbsize 500000000 '+
                    '-use_sw_tback -max_target_seqs 1 -matrix %s'%self.matrix)
        
        if not self.useDisk:
            self.blaster.query = self.query
            self.queryreturn = ''
        else:
            self.queryreturn = self.query + '_' + str(self.uniqueid) + '_return'
    
    def _firstRun(self):
        if self.short:
            res = self.blaster.runBlast(self.query, self.target, self.out,
                         evalue = self.evalue,
                         task='blastp-short',
                         additional=self.additional)
        else:
            res = self.blaster.runBlast(self.query, self.target, self.out,
                         evalue = self.evalue,
                         additional = self.additional)
        
        return res
    
    def _secondRun(self, hit_len = None):
        # Second Blast run
        if not hit_len:
            if self.short:
                hit_len = 29
            else:
                hit_len = 100
                
        if not self.useDisk:
            self.blaster.query = self.blaster.retrieved
                
        if hit_len < 30:
            res = self.blaster.runBlast(self.queryreturn, self.source, self.out,
                     evalue = self.evalue,
                     task='blastp-short',
                     additional=self.additional)
        else:
            res = self.blaster.runBlast(self.queryreturn, self.source, self.out,
                     evalue = self.evalue,
                     additional=self.additional)
            
        return res
    
    def __call__(self):
        if not self.kegg:
            # First Blast run
            res = self._firstRun()
            
            if not res:
                if self.useDisk:
                    try:
                        os.remove(self.out)
                    except:pass
                    
                return (None, self.targetorg, False)
            
            self.blaster.parseBlast(self.out)
            for hits in self.blaster.getHits(self.evalue):
                if len(hits) == 0:
                    break
                targethit = hits[0]
    
                if not self.blaster.retrieveFromDB(self.target, targethit.hit,
                                      out=self.queryreturn):
    
                    if self.useDisk:
                        try:
                            os.remove(self.out)
                        except:pass
                        
                    return (None, self.targetorg, False)
    
                # Second Blast run            
                res = self._secondRun(targethit.hit_len)
                break
        else:
            if not self.blaster.retrieveFromDB(self.target, self.ko_entry,
                                      out=self.queryreturn):
                if self.useDisk:
                    try:
                        os.remove(self.out)
                    except:
                        pass
                return (None, self.targetorg, False)
            res = self._secondRun()
        
        if not res:
        
            if self.useDisk:
                try:
                    os.remove(self.out)
                    os.remove(self.queryreturn)
                except:pass
            return (None, self.targetorg, False)
        
        self.blaster.parseBlast(self.out)
        for hits in self.blaster.getHits(self.evalue):
            if len(hits) == 0:
                return (None, self.targetorg, True)
            sourcehit = hits[0]
            if self.queryid == sourcehit.hit:
                if self.useDisk:
                    os.remove(self.out)
                    os.remove(self.queryreturn)
                if self.kegg:
                    return (self.ko_id,self.queryid, True)
                else:
                    return (sourcehit.query_id.replace('lcl|',''),
                        self.targetorg, True)
            else:
                if self.useDisk:
                    os.remove(self.out)
                    os.remove(self.queryreturn)
                return (None, self.targetorg, True)

        if self.useDisk:
            os.remove(self.out)
        return (None, self.targetorg, True)
