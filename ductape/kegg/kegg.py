#!/usr/bin/env python
"""
Kegg

Kegg Library

Kegg data fetching
"""
import urllib2 as urllib
from ductape.common.commonthread import CommonThread
from ductape.common.utils import get_span
from ductape.common.utils import isOnline
from ductape.kegg.web import kheader
import Queue
import logging
import os
import shutil
import threading
import time
import random

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('ductape.kegg')

################################################################################
# Constants

avoidedPaths = set(['path:rn01110','path:rn01100','path:rn01120',
                'path:ko01100','path:ko01110','path:ko01120',
                'path:map01110','path:map01100','path:map01120',
                'rn01110','rn01100','rn01120',
                'ko01100','ko01110','ko01120',
                'map01110','map01100','map01120',])

################################################################################
# Classes

class MapParser(object):
    '''
    Takes an HTML of a Kegg map and returns a list of ready-to-use html lines
    '''
    def __init__(self,content):
        self.html = content
        self.map = []
        self.parse()
        
    def parse(self):
        b = False
        for l in self.html.split('\n'):
            l = l.lstrip().rstrip().strip()
            if '<map' in l[:4]:
                b = True
                self.map.append(l)
            elif '<area' in l[:6] and b:
                self.map.append(l)
            elif '</map' in l[:5] and b:
                self.map.append(l)
                break
            
        return self.map

class KeggAPI(object):
    '''
    Class KeggAPI
    Connects to KEGG API and performs various tasks
    Fail-safe: if a request fails it tries again and again
    All the results are stored in the attribute result, as well as the inputs
    
    Docs:
    http://www.kegg.jp/kegg/docs/keggapi.html
    http://www.kegg.jp/kegg/rest/weblink.html
    '''
    def __init__(self, keeptrying=False):
        self.baseurl = 'http://www.kegg.jp/'
        self._apiurl = 'http://rest.kegg.jp/'
        self._maplink = 'http://www.kegg.jp/kegg-bin/show_pathway?'
        
        self.failed = False
        
        self.keeptrying = keeptrying
        
        self.clean()
    
    def clean(self):
        self.input = None
        self.result = None
        self.failed = False

    def getEntryTag(self, entry, tag):
        '''
        Get the tag content inside a kegg entry (flat file)
        '''
        b = False
        res = ''
        for line in entry.split('\n'):
            if line.startswith(tag):
                res += line.rstrip().lstrip(tag).lstrip()
                b = True
            elif b and line[0] == ' ':
                res += ' ' + line.rstrip().lstrip()
            elif b and line[0] != ' ':
                b = False
                return res
        
        if res == '':
            return None
        return res.lstrip()
    
    def getLinkTag(self, entry, tag):
        '''
        Get the tag content inside a kegg entry (flat file)
        This variant function extract links from an entry
        '''
        b = False
        for line in entry.split('\n'):
            if line.startswith(tag):
                yield line.rstrip().lstrip(tag).lstrip().split()[-1]
                b = True
            elif b and line[0] == ' ':
                yield line.rstrip().lstrip().split()[-1]
            elif b and line[0] != ' ':
                b = False
        
    def parseLinks(self, links):
        '''
        Parse the results of 
        '''
        d = {}
        for line in links.split('\n'):
            if line == '' or '\t' not in line:continue
            
            k, v = line.split('\t')
            
            if k not in d:
                d[k] = []
            d[k].append(v)
        
        if len(d) == 0:
            return None
        return d

    def getRelease(self, dbver):
        '''
        Takes a string like "Release 64.0+/11-13, Nov 12"
        and returns the release number, as a float (useful for comparison)
        if it fails returns None
        '''
        try:
            s = dbver.split(' ')
            release = s[1]

            while True:
                try:
                    release = float(release)
                    return release
                except:
                    release = release[:-1]
        except:
            logger.debug('Could not parse KEGG database version (%s)'%dbver)
            return None

    def getDBVersion(self, retries=8):
        '''
        Get the KEGG DB version
        Returns a tuple: full version string , release number (None if unparsable)
        '''
        attempts = 0
        while True:
            try:
                self.input = None
                logger.debug('Looking for KEGG db version')
                url = self._apiurl + urllib.quote('info/kegg')
                data = urllib.urlopen(url, timeout=20).read().split('\n')
                
                line = data[1].split('             ')[1]
                self.result = (line, self.getRelease(line))

                return
            except Exception as e:
                attempts += 1
                logger.debug('info failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep((2 + random.random())*attempts)
                try:
                    logger.debug(url)
                except:pass
                if self.keeptrying:continue
                if attempts >= retries:
                    self.failed = True
                    logger.warning('info failed!')
                    return
    
    def getTitle(self, entries, otherTags=[], retries=8):
        '''
        Get the title of a specific KEGG object
        Default behaviour is to return NAME and DEFINITION tags
        Additional tags can be provided in otherTags
        '''
        attempts = 0
        while True:
            try:
                self.input = entries
                logger.debug('Looking for title for %d KEGG entries'%len(entries))
                url = ''
                for entry in entries:
                    url += '%s+'%entry
                    
                # Dummy entry to avoid a rare bug when all the provided entries
                if 'cpd:C00099' not in entries or 'C00099' not in entries:
                    url += 'cpd:C00099'
                #
                
                url = url.rstrip('+')
                url = self._apiurl + 'get/' + urllib.quote(url)
                data = urllib.urlopen(url, timeout=20).read()
                
                self.result = {}
                for lines in data.split('///'):
                    if len(lines) == 1:continue
                    try:
                        shortID = self.getEntryTag(lines,'ENTRY').split(' ')[0]
                        for longID in self.input:
                            if shortID in longID:
                                self.result[longID] = [
                                            self.getEntryTag(lines, 'NAME'),
                                            self.getEntryTag(lines, 'DEFINITION')
                                            ]

                                # TODO: a more general approach for this                                             
                                for tag in otherTags:
                                    value = self.getEntryTag(lines, tag)
                                    self.result[longID].append(value)
                    except:
                        continue
                    
                # Check that every input has a result
                for entry in entries:
                    if entry not in self.result:
                        self.result[entry] = ['','']
                
                return
            except Exception as e:
                attempts += 1
                logger.debug('get failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep((2 + random.random())*attempts)
                try:
                    logger.debug(url)
                except:pass
                if self.keeptrying:continue
                if attempts >= retries:
                    self.failed = True
                    logger.warning('get failed!')
                    return
          
    def getRPair(self, entries, retries=8):
        '''
        Similar to getTitle, but targeting specific features of RPair
        Specifically the two co_id involved and the kind of interaction
        '''
        attempts = 0
        while True:
            try:
                self.input = entries
                logger.debug('Looking for details on %d RPair entries'%len(entries))
                url = ''
                for entry in entries:
                    url += '%s+'%entry
                
                # Dummy entry to avoid a rare bug when all the provided entries
                if 'rp:RP00001' not in entries or 'RP00001' not in entries:
                    url += 'rp:RP00001'
                #
                    
                url = url.rstrip('+')
                url = self._apiurl + 'get/' + urllib.quote(url)
                data = urllib.urlopen(url, timeout=20).read()
                
                self.result = {}
                for lines in data.split('///'):
                    if len(lines) == 1:continue
                    try:
                        shortID = self.getEntryTag(lines,'ENTRY').split(' ')[0]
                        for longID in self.input:
                            if shortID in longID:
                                co1, co2 = self.getEntryTag(lines, 'NAME').split('_')
                                #
                                if not co1.startswith('cpd:'):
                                    co1 = 'cpd:' + co1
                                if not co2.startswith('cpd:'):
                                    co2 = 'cpd:' + co2
                                #
                                kind = self.getEntryTag(lines, 'TYPE')
                                self.result[longID] = [co1,co2,kind]
                    except:
                        continue
                    
                # Check that every input has a result
                for entry in entries:
                    if entry not in self.result:
                        self.result[entry] = ['','','']
                 
                return
            except Exception as e:
                attempts += 1
                logger.debug('get (rpair) failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep((2 + random.random())*attempts)
                try:
                    logger.debug(url)
                except:pass
                if self.keeptrying:continue
                if attempts >= retries:
                    self.failed = True
                    logger.warning('get (rpair) failed!')
                    return
    
    def getIDListFromDB(self, db='pathway', retries=8):
        '''
        Get all the IDs from a specific database
        
        Default: pathway
        '''
        attempts = 0
        while True:
            try:
                self.input = db
                logger.debug('Looking for KEGG IDs from db %s'%db)
                url = self._apiurl + 'list/%s/'%urllib.quote(db)
                
                data = urllib.urlopen(url, timeout=20).read()
                self.result = set([x.split('\t')[0] for x in data.split('\n')])
                try:
                    self.result.remove('')
                except:pass
                self.result = list(self.result)
                return
            except Exception as e:
                attempts += 1
                logger.debug('list (%s) failed! Attempt %d'
                              %(db,attempts))
                logger.debug('%s'%str(e))
                time.sleep((2 + random.random())*attempts)
                try:
                    logger.debug(url)
                except:pass
                if self.keeptrying:continue
                if attempts >= retries:
                    self.failed = True
                    logger.warning('list (%s) failed!'%db)
                    return
    
    def getReactions(self, ko_ids, retries=8):
        '''
        Get the reaction IDs for a given KO list
        '''
        attempts = 0
        while True:
            try:
                self.input = ko_ids
                logger.debug('Looking for KEGG reactions from %d KO IDs'%len(ko_ids))
                url = ''
                for ko_id in ko_ids:
                    url += '%s+'%ko_id
                url = url.rstrip('+')
                
                url = self._apiurl + 'link/reaction/' + urllib.quote(url)
                
                data = urllib.urlopen(url, timeout=20).read()
                self.result = self.parseLinks(data)
                return
            except Exception as e:
                attempts += 1
                logger.debug('link (reaction) failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep((2 + random.random())*attempts)
                try:
                    logger.debug(url)
                except:pass
                if self.keeptrying:continue
                if attempts >= retries:
                    self.failed = True
                    logger.warning('link (reaction) failed!')
                    return
                
    def getPathways(self, re_ids, retries=8):
        '''
        Get the pathway IDs for a given reaction list
        '''
        attempts = 0
        while True:
            try:
                self.input = re_ids
                logger.debug('Looking for KEGG pathways from %d RE IDs'%len(re_ids))
                url = ''
                for re_id in re_ids:
                    url += '%s+'%re_id
                url = url.rstrip('+')
                
                url = self._apiurl + 'link/pathway/' + urllib.quote(url)
                
                data = urllib.urlopen(url, timeout=20).read()
                self.result = self.parseLinks(data)
                return
            except Exception as e:
                attempts += 1
                logger.debug('link (pathway) failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep((2 + random.random())*attempts)
                try:
                    logger.debug(url)
                except:pass
                if self.keeptrying:continue
                if attempts >= retries:
                    self.failed = True
                    logger.warning('link (pathway) failed!')
                    return
    
    def getReactionsByComp(self, co_ids, retries=8):
        '''
        Get the reactions IDs for a given compound list
        '''
        attempts = 0
        while True:
            try:
                self.input = co_ids
                logger.debug('Looking for KEGG reactions from %d CO IDs'%len(co_ids))
                url = ''
                for co_id in co_ids:
                    url += '%s+'%co_id
                url = url.rstrip('+')
                
                url = self._apiurl + 'link/reaction/' + urllib.quote(url)
                
                data = urllib.urlopen(url, timeout=20).read()
                self.result = self.parseLinks(data)
                return
            except Exception as e:
                attempts += 1
                logger.debug('link (reaction) failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep((2 + random.random())*attempts)
                try:
                    logger.debug(url)
                except:pass
                if self.keeptrying:continue
                if attempts >= retries:
                    self.failed = True
                    logger.warning('link (reaction) failed!')
                    return
    
    def getReactionsFromPath(self, path_ids, retries=8):
        '''
        Get the reaction IDs for a given pathway list
        '''
        attempts = 0
        while True:
            try:
                self.input = path_ids
                logger.debug('Looking for KEGG reactions from %d PATH IDs'%len(path_ids))
                url = ''
                for path_id in path_ids:
                    url += '%s+'%path_id
                url = url.rstrip('+')
                
                url = self._apiurl + 'link/reaction/' + urllib.quote(url)
                
                data = urllib.urlopen(url, timeout=20).read()
                self.result = self.parseLinks(data)
                return
            except Exception as e:
                attempts += 1
                logger.debug('link (reaction) failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep((2 + random.random())*attempts)
                try:
                    logger.debug(url)
                except:pass
                if self.keeptrying:continue
                if attempts >= retries:
                    self.failed = True
                    logger.warning('link (reaction) failed!')
                    return
    
    def getRPairsFromReaction(self, entries, retries=8):
        '''
        Get the rpair IDs for a given reaction list
        '''
        attempts = 0
        while True:
            try:
                self.input = entries
                logger.debug('Looking for RClass for %d KEGG entries'%len(entries))
                url = ''
                for entry in entries:
                    url += '%s+'%entry
                    
                # Dummy entry to avoid a rare bug when all the provided entries
                if 'cpd:C00099' not in entries or 'C00099' not in entries:
                    url += 'cpd:C00099'
                #
                
                url = url.rstrip('+')
                url = self._apiurl + 'get/' + urllib.quote(url)
                data = urllib.urlopen(url, timeout=20).read()
                
                self.result = {}
                for lines in data.split('///'):
                    if len(lines) == 1:continue
                    try:
                        shortID = self.getEntryTag(lines,'ENTRY').split(' ')[0]
                        for longID in self.input:
                            if shortID in longID:
                                for rclass in self.getLinkTag(lines, 'RCLASS'):
                                    self.result[longID] = self.result.get(longID,
                                                                          set())
                                    self.result[longID].add(rclass)
                    except:
                        continue
                    
                return
            except Exception as e:
                attempts += 1
                logger.debug('link (rpair) failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep((2 + random.random())*attempts)
                try:
                    logger.debug(url)
                except:pass
                if self.keeptrying:continue
                if attempts >= retries:
                    self.failed = True
                    logger.warning('link (rpair) failed!')
                    return
                
    def getCompoundsFromReaction(self, re_ids, retries=8):
        '''
        Get the compound IDs for a given reaction list
        '''
        attempts = 0
        while True:
            try:
                self.input = re_ids
                logger.debug('Looking for KEGG compounds from %d RE IDs'%len(re_ids))
                url = ''
                for re_id in re_ids:
                    url += '%s+'%re_id
                url = url.rstrip('+')
                
                url = self._apiurl + 'link/compound/' + urllib.quote(url)
                
                data = urllib.urlopen(url, timeout=20).read()
                self.result = self.parseLinks(data)
                return
            except Exception as e:
                attempts += 1
                logger.debug('link (compound) failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep((2 + random.random())*attempts)
                try:
                    logger.debug(url)
                except:pass
                if self.keeptrying:continue
                if attempts >= retries:
                    self.failed = True
                    logger.warning('link (compound) failed!')
                    return
    
    def getCompoundsFromPath(self, path_ids, retries=8):
        '''
        Get the compound IDs for a given pathway list
        '''
        attempts = 0
        while True:
            try:
                self.input = path_ids
                logger.debug('Looking for KEGG compounds from %d PATH IDs'%len(path_ids))
                url = ''
                for path_id in path_ids:
                    url += '%s+'%path_id
                url = url.rstrip('+')
                
                url = self._apiurl + 'link/compound/' + urllib.quote(url)
                
                data = urllib.urlopen(url, timeout=20).read()
                self.result = self.parseLinks(data)
                return
            except Exception as e:
                attempts += 1
                logger.debug('link (compound) failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep((2 + random.random())*attempts)
                try:
                    logger.debug(url)
                except:pass
                if self.keeptrying:continue
                if attempts >= retries:
                    self.failed = True
                    logger.warning('link (compound) failed!')
                    return
    
    def getHTMLColoredPathway(self, path_id, obj_list, color_list,
                                    border_list=None, retries=8):
        '''
        Get the URL of the colored pathway and return its content
        If it fails, an exception is thrown
        '''
        attempts = 0
        while True:
            try:
                self.input = path_id
                
                # Fix color codes
                for i in range(len(color_list)):
                    if '#' in color_list[i]:
                        color_list[i] = color_list[i].replace('#', '%23')
                #
                
                logger.debug('Looking for KEGG colored map from %s'%path_id)
                url = path_id.lstrip('path:') + '/default%3white/'
                for i in range(len(obj_list)):
                    if border_list is not None and border_list[i] is not None:
                        url += obj_list[i] + '%09' + color_list[i] + ',' + border_list[i] + '/'
                    else:
                        url += obj_list[i] + '%09' + color_list[i] + '/'
               
                # Cannot quote this url, no colored pathway can be obtained then
                #url = self._maplink + urllib.quote(url)
                url = self._maplink + url
                
                logger.debug(url)
                
                sock=urllib.urlopen(url, timeout=60)
                self.result = sock.read()
                sock.close()
                return
            except Exception as e:
                attempts += 1
                logger.debug('show_pathway failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep((2 + random.random())*attempts)
                if self.keeptrying:continue
                if attempts >= retries:
                    self.failed = True
                    logger.warning('show_pathway failed!')
                    return

class KeggColor(object):
    '''
    Class KeggColor
    Holds the color information to be passed to MapsFetcher
    One object for each pathway
    '''
    def __init__(self, path, htmlmap= '', reactions={}, compounds={},
                 borders={}):
        self.path = path 
        self.htmlmap = htmlmap
        self.reactions = reactions
        self.compounds = compounds
        # Objects that need to have a coloured border
        self.borders = borders
    
    def setMap(self, htmlmap):
        self.htmlmap = htmlmap
    
    def setReactions(self, reactions):
        self.reactions = reactions
        
    def setCompounds(self, compounds):
        self.compounds = compounds
        
    def setBorders(self, borders):
        self.borders = borders
        
    def getAll(self):
        '''
        Returns a tuple --> objects, color
        '''
        objs = [x for x in self.reactions]
        colors = [self.reactions[x] for x in self.reactions]
        objs += [x for x in self.compounds]
        colors += [self.compounds[x] for x in self.compounds]
                
        return objs,colors
    
    def getBorders(self):
        '''
        Returns a tuple --> objects, color
        '''
        objs = [x for x in self.reactions]
        objs += [x for x in self.compounds]
        
        colors = []
        for x in self.reactions.keys()+self.compounds.keys():
            if x in self.borders.keys():
                colors.append(self.borders[x])
            else:
                colors.append(None)
                
        return objs,colors
        
class KeggDetails(object):
    '''
    Class KoDetails
    All the informations returned by Mappers are contained here
    '''
    def __init__(self):
        # Details
        self.ko = None
        self.react = None
        self.comp = None
        self.path = None
        self.rpair = None
        # Links
        self.koreact = None
        self.pathreact = None
        self.pathcomp = None
        self.compreact = None
        self.reactcomp = None
        self.reactrpair = None
        self.rpairreact = None
        # Maps
        self.pathmaps = None
    
    def _purgeDetails(self,det):
        erase = []
        if not det:
            return det
        
        for key, value in det.iteritems():
            if not value:
                erase.append(key)
        
        for key in erase:
            del det[key]
        
        return det
    
    def setDetails(self, ko=None, react=None, comp=None, path=None, rpair=None):
        self.ko = self._purgeDetails(ko)
        self.react = self._purgeDetails(react)
        self.comp = self._purgeDetails(comp)
        self.path = self._purgeDetails(path)
        self.rpair = self._purgeDetails(rpair)
        
    def setLinks(self, koreact=None, pathreact=None, pathcomp=None,
                 compreact=None, reactcomp=None, reactrpair=None,
                 rpairreact=None):
        self.koreact = {}
        if koreact:
            for k,v in koreact.iteritems():
                self.koreact[k] = []
                for i in v:
                    self.koreact[k].append(str(i))

        self.pathreact = {}
        if pathreact:
            for k,v in pathreact.iteritems():
                self.pathreact[k] = []
                for i in v:
                    self.pathreact[k].append(str(i))
                    
        self.pathcomp = {}
        if pathcomp:
            for k,v in pathcomp.iteritems():
                self.pathcomp[k] = []
                for i in v:
                    self.pathcomp[k].append(str(i))
                    
        self.compreact = {}
        if compreact:
            for k,v in compreact.iteritems():
                self.compreact[k] = []
                for i in v:
                    self.compreact[k].append(str(i))
                    
                self.compreact = {}
        
        self.reactcomp = {}
        if reactcomp:
            for k,v in reactcomp.iteritems():
                self.reactcomp[k] = []
                for i in v:
                    self.reactcomp[k].append(str(i))
                    
        self.reactrpair = {}
        if reactrpair:
            for k,v in reactrpair.iteritems():
                self.reactrpair[k] = []
                for i in v:
                    self.reactrpair[k].append(str(i))
                    
        self.rpairreact = {}
        if rpairreact:
            for k,v in rpairreact.iteritems():
                self.rpairreact[k] = []
                for i in v:
                    self.rpairreact[k].append(str(i))
    
    def setMaps(self, maps):
        self.pathmaps = maps
    
    def getKO(self):
        return self.ko
    
    def getReact(self):
        return self.react
    
    def getComp(self):
        return self.comp
    
    def getPath(self):
        return self.path
    
    def getKOLinks(self):
        return self.koreact
    
    def getCompLinks(self):
        return self.comppath
    
    def getPathLinks(self):
        return self.pathreact, self.pathcomp
    
    def getMaps(self):
        return self.pathmaps

class BaseKegg(CommonThread):
    def __init__(self, threads=40, keeptrying=False, queue=Queue.Queue()):
        CommonThread.__init__(self,queue)
        
        # Kegg connection
        self.handlers = []
        self._hindex = 0
        self.numThreads = threads
        
        for i in range(self.numThreads):
            obj = KeggAPI(keeptrying)
            self.handlers.append(obj)
            
        self.cleanHandlers()
    
    def cleanHandlers(self):
        for handler in self.handlers:
            handler.clean()
        self._hindex = 0
    
    def getHandler(self):
        if self._hindex >= len(self.handlers):
            self._hindex = 0
        handler = self.handlers[self._hindex] 
        self._hindex += 1
        return handler
        
    def checkConnection(self):
        '''
        Check if there are connection problems
        First check the two IP addresses, then the URL
        '''
        check = [KeggAPI().baseurl, KeggAPI()._apiurl]
        online = False
        for addr in check:
            try:
                isOnline(addr)
                online = True
            except:
                logger.debug('address %s not working'%addr)        
        if not online:
            raise Exception('KEGG seems to be offline')
            
class BaseMapper(BaseKegg):
    def __init__(self, threads=40, avoid=[], keeptrying=False,
                        queue=Queue.Queue()):
        BaseKegg.__init__(self, threads=threads, keeptrying=keeptrying,
                                queue=queue)

        # Skip these
        self.avoid = avoid
        
        # Results
        self.reactdet = {}
        self.rpairdet = {}
        self.pathdet = {}
        self.pathreact = {}
        self.pathcomp = {}
        self.pathmap = {}
        self.compdet = {}
        self.reactpath = {}
        self.reactcomp = {}
        self.compreact = {}
        self.rpairreact = {}
        self.reactrpair = {}
        
        # Output
        self.result = None
        
    def getReactDetails(self):
        pieces = [p for p in get_span(self.reactdet.keys(), 9)]
        for piece in get_span(pieces, self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += len([i for p in piece for i in p])
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for ids in piece:
                remove = set()
                for i in ids:
                    if i in self.avoid:
                        remove.add(i)
                for i in remove:
                    ids.remove(i)
                
                if len(ids) == 0:
                    continue
                
                obj = threading.Thread(
                                target = self.getHandler().getTitle,
                                args = (ids,['ENZYME'],))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            if len(threads) == 0:
                continue
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if handler.failed:
                    logger.error('KEGG API error, aborting')
                    raise IOError('KEGG API error')
                
                if not handler.result:
                    logger.debug('Found an empty handler')
                    continue
                
                for kid, title in handler.result.iteritems():
                    self.reactdet[kid] = title
    
    def getRPairDetails(self):
        pieces = [p for p in get_span(self.rpairdet.keys(), 9)]
        for piece in get_span(pieces, self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += len([i for p in piece for i in p])
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            for ids in piece:
                remove = set()
                for i in ids:
                    if i in self.avoid:
                        remove.add(i)
                for i in remove:
                    ids.remove(i)
                    
                if len(ids) == 0:
                    continue
                
                for rid in ids:
                    self.rpairdet[rid] = [rid.split('_')[0],
                                          rid.split('_')[1],
                                          'main']
    
    def getPathDetails(self):
        pieces = [p for p in get_span(self.pathdet.keys(), 9)]
        for piece in get_span(pieces, self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += len([i for p in piece for i in p])
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for ids in piece:
                remove = set()
                for i in ids:
                    if i in self.avoid:
                        remove.add(i)
                for i in remove:
                    ids.remove(i)
                    
                if len(ids) == 0:
                    continue
                
                obj = threading.Thread(
                                target = self.getHandler().getTitle,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            if len(threads) == 0:
                continue
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if handler.failed:
                    logger.error('KEGG API error, aborting')
                    raise IOError('KEGG API error')
                
                if not handler.result:
                    logger.debug('Found an empty handler')
                    continue
                
                for kid, title in handler.result.iteritems():
                    self.pathdet[kid] = title
    
    def getMapsDetails(self):
        for piece in get_span(self.pathdet.keys(), self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += self.numThreads
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for path in piece:
                if path in self.avoid:
                    continue
                
                obj = threading.Thread(
                            target = self.handlers[piece.index(path)].getHTMLColoredPathway,
                        args = (path,[],[],))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            if len(threads) == 0:
                continue
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if handler.failed:
                    logger.error('KEGG API error, aborting')
                    raise IOError('KEGG API error')
                
                if not handler.result:
                    logger.debug('Found an empty handler')
                    continue
                parser = MapParser(handler.result)
                self.pathmap[handler.input] = parser.map
    
    def getPathReactions(self):
        pieces = [p for p in get_span(self.pathdet.keys(), 80)]
        for piece in get_span(pieces, self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += len([i for p in piece for i in p])
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for ids in piece:
                remove = set()
                for i in ids:
                    if i in self.avoid:
                        remove.add(i)
                for i in remove:
                    ids.remove(i)
                    
                if len(ids) == 0:
                    continue
                
                obj = threading.Thread(
                                target = self.getHandler().getReactionsFromPath,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            if len(threads) == 0:
                continue
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if handler.failed:
                    logger.error('KEGG API error, aborting')
                    raise IOError('KEGG API error')
                
                if not handler.result:
                    logger.debug('Found an empty handler')
                    continue
                
                for path, reacts in handler.result.iteritems():
                    if path not in self.pathreact:
                        self.pathreact[path] = reacts
                reacts = set([v for vs in handler.result.itervalues() for v in vs])
                for react in reacts:
                    if react not in self.reactdet:
                        self.reactdet[react] = None
                        
    def getPathCompounds(self):
        pieces = [p for p in get_span(self.pathdet.keys(), 80)]
        for piece in get_span(pieces, self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += len([i for p in piece for i in p])
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for ids in piece:
                remove = set()
                for i in ids:
                    if i in self.avoid:
                        remove.add(i)
                for i in remove:
                    ids.remove(i)
                    
                if len(ids) == 0:
                    continue
                
                obj = threading.Thread(
                                target = self.getHandler().getCompoundsFromPath,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            if len(threads) == 0:
                continue
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if handler.failed:
                    logger.error('KEGG API error, aborting')
                    raise IOError('KEGG API error')
                
                if not handler.result:
                    logger.debug('Found an empty handler')
                    continue
                
                for path, comps in handler.result.iteritems():
                    if path not in self.pathcomp:
                        self.pathcomp[path] = comps
                comps = set([v for vs in handler.result.itervalues() for v in vs])
                for comp in comps:
                    if comp not in self.compdet:
                        self.compdet[comp] = None
                        
    def getCompDetails(self):
        pieces = [p for p in get_span(self.compdet.keys(), 9)]
        for piece in get_span(pieces, self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += len([i for p in piece for i in p])
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for ids in piece:
                remove = set()
                for i in ids:
                    if i in self.avoid:
                        remove.add(i)
                for i in remove:
                    ids.remove(i)
                    
                if len(ids) == 0:
                    continue
                
                obj = threading.Thread(
                                target = self.getHandler().getTitle,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            if len(threads) == 0:
                continue
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if handler.failed:
                    logger.error('KEGG API error, aborting')
                    raise IOError('KEGG API error')
                
                if not handler.result:
                    logger.debug('Found an empty handler')
                    continue
                
                for kid, title in handler.result.iteritems():
                    self.compdet[kid] = title
    
    def getPathways(self):
        pieces = [p for p in get_span(self.reactdet.keys(), 80)]
        for piece in get_span(pieces, self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += len([i for p in piece for i in p])
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for ids in piece:
                remove = set()
                for i in ids:
                    if i in self.avoid:
                        remove.add(i)
                for i in remove:
                    ids.remove(i)
                    
                if len(ids) == 0:
                    continue
                
                obj = threading.Thread(
                                target = self.getHandler().getPathways,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            if len(threads) == 0:
                continue
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if handler.failed:
                    logger.error('KEGG API error, aborting')
                    raise IOError('KEGG API error')
                
                if not handler.result:
                    logger.debug('Found an empty handler')
                    continue
                
                for react, paths in handler.result.iteritems():
                    if react not in self.reactpath:
                        self.reactpath[react] = []
                    for path in paths:
                        if path.startswith('path:map'):continue
                        self.reactpath[react].append(path)
                paths = set([v for vs in handler.result.itervalues() for v in vs])
                for path in paths:
                    if path not in self.pathdet and not path.startswith('path:map'):
                        self.pathdet[path] = None
                        
    def getReactCompounds(self):
        pieces = [p for p in get_span(self.reactdet.keys(), 80)]
        for piece in get_span(pieces, self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += len([i for p in piece for i in p])
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for ids in piece:
                remove = set()
                for i in ids:
                    if i in self.avoid:
                        remove.add(i)
                for i in remove:
                    ids.remove(i)
                    
                if len(ids) == 0:
                    continue
                
                obj = threading.Thread(
                                target = self.getHandler().getCompoundsFromReaction,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            if len(threads) == 0:
                continue
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if handler.failed:
                    logger.error('KEGG API error, aborting')
                    raise IOError('KEGG API error')
                
                if not handler.result:
                    logger.debug('Found an empty handler')
                    continue
                
                for react, comps in handler.result.iteritems():
                    if react not in self.reactcomp:
                        self.reactcomp[react] = comps
                comps = set([v for vs in handler.result.itervalues() for v in vs])
                for comp in comps:
                    if comp not in self.compdet:
                        self.compdet[comp] = None
                        
    def getCompoundReacts(self):
        pieces = [p for p in get_span(self.compdet.keys(), 80)]
        for piece in get_span(pieces, self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += len([i for p in piece for i in p])
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for ids in piece:
                remove = set()
                for i in ids:
                    if i in self.avoid:
                        remove.add(i)
                for i in remove:
                    ids.remove(i)
                    
                if len(ids) == 0:
                    continue
                
                obj = threading.Thread(
                                target = self.getHandler().getReactionsByComp,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            if len(threads) == 0:
                continue
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if handler.failed:
                    logger.error('KEGG API error, aborting')
                    raise IOError('KEGG API error')
                
                if not handler.result:
                    logger.debug('Found an empty handler')
                    continue
                
                for comp, reacts in handler.result.iteritems():
                    if comp not in self.compreact:
                        self.compreact[comp] = reacts
                reacts = set([v for vs in handler.result.itervalues() for v in vs])
                for react in reacts:
                    if react not in self.reactdet:
                        self.reactdet[react] = None
    
    def getReactRPairs(self):
        pieces = [p for p in get_span(self.reactdet.keys(), 80)]
        for piece in get_span(pieces, self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += len([i for p in piece for i in p])
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for ids in piece:
                remove = set()
                for i in ids:
                    if i in self.avoid:
                        remove.add(i)
                for i in remove:
                    ids.remove(i)
                    
                if len(ids) == 0:
                    continue
                
                obj = threading.Thread(
                                target = self.getHandler().getRPairsFromReaction,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            if len(threads) == 0:
                continue
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if handler.failed:
                    logger.error('KEGG API error, aborting')
                    raise IOError('KEGG API error')
                
                if not handler.result:
                    logger.debug('Found an empty handler')
                    continue
                
                for react, rpairs in handler.result.iteritems():
                    if react not in self.reactrpair:
                        self.reactrpair[react] = rpairs
                rpairs = set([v for vs in handler.result.itervalues() for v in vs])
                for rpair in rpairs:
                    if rpair not in self.rpairdet:
                        self.rpairdet[rpair] = None

class KoMapper(BaseMapper):
    '''
    From a list of KO IDs returns various details in an object
    KO --> title and details
    KO --> reactions (and titles)
    reactions --> pathways (and titles)
    pathways --> reactions, compounds (with titles)
    maps --> for each pathway, the html maps (!!!)
    '''
    
    _statusDesc = {0:'Not started',
               1:'Checking connectivity',
               2:'Fetching reactions',
               3:'Fetching rpairs',
               4:'Fetching pathways',
               5:'Fetching pathways content',
               6:'Fetching reactions - compounds links',
               7:'Fetching details on KEGG entries',
               8:'Crafting results'}
    
    _substatuses = [2,3,4,5,6,7]
    
    def __init__(self, ko_list, threads=40, avoid=[], keeptrying=False,
                        queue=Queue.Queue()):
        BaseMapper.__init__(self, threads=threads, avoid=avoid,
                            keeptrying=keeptrying, queue=queue)
        # Kegg
        self.ko = ko_list
        
        # Results
        self.kodet = {}
        self.koreact = {}
    
    def getKOdet(self):
        pieces = [p for p in get_span(self.ko, 9)]
        for piece in get_span(pieces, self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += len([i for p in piece for i in p])
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for ids in piece:
                remove = set()
                for i in ids:
                    if i in self.avoid:
                        remove.add(i)
                for i in remove:
                    ids.remove(i)
                    
                if len(ids) == 0:
                    continue
                
                obj = threading.Thread(
                                target = self.getHandler().getTitle,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            if len(threads) == 0:
                continue
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if handler.failed:
                    logger.error('KEGG API error, aborting')
                    raise IOError('KEGG API error')
                
                if not handler.result:
                    logger.debug('Found an empty handler')
                    continue
                
                for kid, title in handler.result.iteritems():
                    self.kodet[kid] = title
                
    def getReactions(self):
        pieces = [p for p in get_span(self.ko, 80)]
        for piece in get_span(pieces, self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += len([i for p in piece for i in p])
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for ids in piece:
                remove = set()
                for i in ids:
                    if i in self.avoid:
                        remove.add(i)
                for i in remove:
                    ids.remove(i)
                    
                if len(ids) == 0:
                    continue
                
                obj = threading.Thread(
                                target = self.getHandler().getReactions,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            if len(threads) == 0:
                continue
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if handler.failed:
                    logger.error('KEGG API error, aborting')
                    raise IOError('KEGG API error')
                
                if not handler.result:
                    logger.debug('Found an empty handler')
                    continue
                
                for ko, reacts in handler.result.iteritems():
                    if ko not in self.koreact:
                        self.koreact[ko] = reacts
                reacts = set([v for vs in handler.result.itervalues() for v in vs])
                for react in reacts:
                    if react not in self.reactdet:
                        self.reactdet[react] = None
    
    def run(self):
        self.updateStatus()
        try:
            self.checkConnection()
        except Exception as e:
            self.sendFailure(str(e))
            return
    
        # Reactions
        self._maxsubstatus = len(self.ko)
        self.updateStatus()
        try:
            self.getReactions()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Related rpairs
        self._maxsubstatus = len(self.reactdet)
        self.updateStatus()
        logger.warning('Using RCLASS attribute of KEGG reactions, as the RPAIR database has now been discontinued')
        try:
            self.getReactRPairs()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Related pathways...
        self._maxsubstatus = len(self.reactdet)
        self.updateStatus()
        try:
            self.getPathways()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Pathways contents...
        # 1. Reactions
        self._maxsubstatus = len(self.pathdet)
        self.updateStatus()
        try:
            self.getPathReactions()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # 2. Compounds
        self._maxsubstatus = len(self.pathdet)
        try:
            self.getPathCompounds()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Compounds for each reaction
        self._maxsubstatus = len(self.reactdet)
        self.updateStatus()
        try:
            self.getReactCompounds()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Reactions for each compound
        self._maxsubstatus = len(self.compdet)
        try:
            self.getCompoundReacts()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # KO details
        self._maxsubstatus = len(self.ko)
        self.updateStatus()
        try:
            self.getKOdet()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Pathway details
        self._maxsubstatus = len(self.pathdet)
        try:
            self.getPathDetails()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Pathway HTML maps (!!!)
        self._maxsubstatus = len(self.pathdet)
        try:
            self.getMapsDetails()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Reaction details
        self._maxsubstatus = len(self.reactdet)
        try:
            self.getReactDetails()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Compound details
        self._maxsubstatus = len(self.compdet)
        try:
            self.getCompDetails()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # RPair details
        self._maxsubstatus = len(self.rpairdet)
        logger.warning('Using RCLASS attribute of KEGG reactions, as the RPAIR database has now been discontinued')
        try:
            self.getRPairDetails()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Prepare the output object
        self.updateStatus()
        self.result = KeggDetails()
        self.result.setDetails(self.kodet, self.reactdet,
                               self.compdet, self.pathdet, self.rpairdet)
        self.result.setLinks(koreact=self.koreact, pathreact=self.pathreact, 
                             pathcomp=self.pathcomp, reactcomp=self.reactcomp,
                             compreact=self.compreact,
                             rpairreact=self.rpairreact,
                             reactrpair=self.reactrpair)
        self.result.setMaps(self.pathmap)

class CompMapper(BaseMapper):
    '''
    From a list of CO IDs returns various details in an object
    CO --> title and details
    CO --> pathways (and titles)
    pathways --> reactions, compounds (with titles)
    maps --> for each pathway, the html maps (!!!)
    '''
    
    _statusDesc = {0:'Not started',
               1:'Checking connectivity',
               2:'Fetching reactions',
               3:'Fetching rpairs',
               4:'Fetching pathways',
               5:'Fetching pathways content',
               6:'Fetching reactions - compounds links',
               7:'Fetching details on KEGG entries',
               8:'Crafting results'}
    
    _substatuses = [2,3,4,5,6,7]
    
    def __init__(self, co_list, threads=40, avoid=[], keeptrying=False,
                        queue=Queue.Queue()):
        BaseMapper.__init__(self, threads=threads, avoid=avoid,
                            keeptrying=keeptrying, queue=queue)
        # Kegg
        self.co = co_list
        
        # Results
        self.comppath = {}
    
    def run(self):
        self.updateStatus()
        try:
            self.checkConnection()
        except Exception as e:
            self.sendFailure(str(e))
            return
    
        # Reactions
        for co_id in self.co:
            self.compdet[co_id] = None
        self._maxsubstatus = len(self.compdet)
        self.updateStatus()
        try:
            self.getCompoundReacts()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Related rpairs
        self._maxsubstatus = len(self.reactdet)
        self.updateStatus()
        logger.warning('Using RCLASS attribute of KEGG reactions, as the RPAIR database has now been discontinued')
        try:
            self.getReactRPairs()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Related pathways...
        self._maxsubstatus = len(self.reactdet)
        self.updateStatus()
        try:
            self.getPathways()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Pathways contents...
        # 1. Reactions
        self._maxsubstatus = len(self.pathdet)
        self.updateStatus()
        try:
            self.getPathReactions()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # 2. Compounds
        self._maxsubstatus = len(self.pathdet)
        try:
            self.getPathCompounds()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Compunds for each reaction
        self._maxsubstatus = len(self.reactdet)
        self.updateStatus()
        try:
            self.getReactCompounds()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Reactions for each compund
        self._maxsubstatus = len(self.compdet)
        try:
            self.getCompoundReacts()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Pathway details
        self._maxsubstatus = len(self.pathdet)
        self.updateStatus()
        try:
            self.getPathDetails()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Pathway HTML maps (!!!)
        self._maxsubstatus = len(self.pathdet)
        try:
            self.getMapsDetails()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Reaction details
        self._maxsubstatus = len(self.reactdet)
        try:
            self.getReactDetails()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Compound details
        self._maxsubstatus = len(self.compdet)
        try:
            self.getCompDetails()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # RPair details
        self._maxsubstatus = len(self.rpairdet)
        logger.warning('Using RCLASS attribute of KEGG reactions, as the RPAIR database has now been discontinued')
        try:
            self.getRPairDetails()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Prepare the output object
        self.updateStatus()
        self.result = KeggDetails()
        self.result.setDetails(react=self.reactdet,
                               comp=self.compdet, path=self.pathdet,
                               rpair=self.rpairdet)
        self.result.setLinks(pathreact=self.pathreact, 
                             pathcomp=self.pathcomp, compreact=self.compreact,
                             reactcomp=self.reactcomp,
                             reactrpair=self.reactrpair,
                             rpairreact=self.rpairreact)
        self.result.setMaps(self.pathmap)

class MapsFetcher(BaseKegg):
    '''
    Class MapsFetcher
    Download colored Kegg maps (png or URLs)
    Input: color_objs (KeggColor list), picture, htmls, urls, prefix, legend (file)
    Output: tuple(list of png filenames, list of HTML files, list of URLs)
    '''
    
    _statusDesc = {0:'Not started',
               1:'Checking connectivity',
               2:'Making room',
               3:'Fetching maps (pictures)',
               4:'Generating interactive web pages'}
    
    _substatuses = [3]
    
    def __init__(self, color_objs, pictures=True, html=True, prefix='', 
                 legend=None, threads=40, keeptrying=False,
                 queue=Queue.Queue()):
        BaseKegg.__init__(self, threads=threads, keeptrying=keeptrying,
                          queue=queue)
        
        self.colors = color_objs
        self.pictures = bool(pictures)
        self.web = bool(html)
        self.legend = legend 
        
        self._keggroom = None
        self._prefix = prefix 
        
        # Outputs
        self.pics = []
        self.webpages = []
        self.pages = []
        self.result = (self.pics, self.webpages, self.pages)
    
    def makeRoom(self,location=''):
        '''
        Creates a tmp directory in the desired location
        '''
              
        # KEGG database path
        try:
            path = os.path.abspath(location)
            path = os.path.join(path, 'tmp')
            try:os.mkdir(path)
            except:pass
            path = os.path.join(path, 'keggmaps')
            try:os.mkdir(path)
            except:pass
            path = os.path.join(path, self._prefix)
            self._keggroom = path
            os.mkdir(path)
        except:
            logger.debug('Temporary directory creation failed! %s'
                          %path)
    
    def copyLegend(self):
        '''Copy the legend in the target directory'''
        if self.legend and os.path.exists(self.legend):
            legend = os.path.join(self._keggroom, 'legend.png')
            shutil.copyfile(self.legend, legend)
            
            return legend
        return None
    
    def getMaps(self):
        legend = self.copyLegend()
        
        for piece in get_span(self.colors, self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += self.numThreads
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for kmap in piece:
                path = kmap.path
                
                # Skip the general maps
                if path in avoidedPaths:
                    logger.debug('Skipping general pathway %s'%path)
                    continue
                #
                
                objs,colors = kmap.getAll()
                dummy,borders = kmap.getBorders()
                
                obj = threading.Thread(
                        target = self.handlers[piece.index(kmap)].getHTMLColoredPathway,
                        args = (path,objs,colors,borders,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            if len(threads) == 0:
                continue
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if handler.failed:
                    logger.error('KEGG API error, aborting')
                    raise IOError('KEGG API error')
                
                if not handler.result:
                    logger.debug('Found an empty handler')
                    continue
                
                fname = os.path.join(self._keggroom,handler.input)
                fname = fname+'.png'
                
                # Fetch the map picture
                # Hoping it won't change much in the future
                for line in handler.result.split('\n'):
                    if ('<img' in line
                        and 'pathwayimage' in line
                        and 'usemap="#mapdata"' in line):
                        urlimage = 'http://www.kegg.jp/' + line.split('src="')[1].split('"')[0]

                        sock=urllib.urlopen(urlimage, timeout=30)
                        pic = sock.read()
                        sock.close()
                        
                        fOut = open(fname,'w')
                        fOut.write(pic)
                        fOut.close()
                        self.pics.append(fname)
                #
    
    def getWebPages(self):
        # TODO: nicer web pages
        legend = self.copyLegend()
        legend = os.path.split(legend)[-1]
        
        if legend:
            fname = os.path.join(self._keggroom,'legend.html')
            
            fOut = open(fname,'w')
            fOut.write('<html>\n<head></head>\n<body>\n')
            fOut.write('''<div align="center">
                        <img src="./%s" />
                        </div>\n'''%
                       (legend))
            fOut.write('</body>\n</html>')
            fOut.close()
        
        for path in self.colors:
            myindex = self.colors.index(path)
            
            logger.debug('Writing interactive web page for %s'%path.path)
            
            fname = os.path.join(self._keggroom,path.path)
            fname = fname+'.html'
            
            fOut = open(fname,'w')
            fOut.write('<html>\n%s\n<body>\n'%kheader)
            
            # Navigation
            fOut.write('''<h2 align="center">\n''')
            fOut.write('''<a href="./%s.html">&laquo;</a>%s'''%(
                       self.colors[myindex-1].path, path.path))
            try:
                next = self.colors[myindex+1]
            except:
                next = self.colors[0]
            fOut.write('''<a href="./%s.html">&raquo;</a>\n</h2>\n'''%
                       (next.path))
            
            if legend:
                fOut.write('''<h3 align="center">
                            <a href="./legend.html">Color scheme</a>
                            </h3>\n''')
            
            fOut.write('''<div align="center">
                        <img src="./%s" usemap="#mapdata" border="0" />
                        </div>\n'''%
                       (path.path+'.png'))
            
            html = path.htmlmap.split('\n')
            newhtml = []
            for line in html:
                line = line.replace('href="/dbget-bin/www_bget?',
               'target="_blank" href="http://www.genome.jp/dbget-bin/www_bget?')
                
                if '/kegg-bin/show_pathway?' in line:
                    s = line.split('/kegg-bin/show_pathway?')
                    s1 = s[1].split('"')
                    s1[0] += '.html'
                    line1 = '"'.join(s1)
                    line = './'.join([s[0]] + [line1])
                
                newhtml.append(line)
                
            fOut.write('%s\n'%'\n'.join(newhtml))
            fOut.write('<div id="poplay" class="poplay" />\n</body>\n</html>')
            fOut.close()
            
            self.webpages.append(fname)
    
    def run(self):
        self.updateStatus()
        try:
            self.checkConnection()
        except Exception as e:
            self.sendFailure(str(e))
            return
    
        self.updateStatus()
        self.makeRoom()
        
        if self.killed:
            return
        
        # ':' bugfix
        # the ':' char causes various problems in windows folders
        for path in self.colors:
            if ':' in path.path:
                path.path = path.path.split(':')[1]
        
        if self.pictures:
            self._maxsubstatus = len(self.colors)
            self.updateStatus()
            try:
                self.getMaps()
            except Exception as e:
                self.sendFailure(e)
                return
            self.cleanHandlers()
            self.resetSubStatus()
        else:
            self.updateStatus(send=False)
            
        if self.killed:
            return
        
        if self.web:
            self.updateStatus()
            try:
                self.getWebPages()
            except Exception as e:
                self.sendFailure(e)
                return
        else:
            self.updateStatus(send=False)

class KeggNet(BaseMapper):
    '''
    Fetch as much details as possible from the KEGG database,
    starting from the pathways list
    '''
    
    _statusDesc = {0:'Not started',
               1:'Checking connectivity',
               2:'Fetching pathways',
               3:'Fetching compounds',
               4:'Fetching compounds - reactions links',
               5:'Fetching reactions',
               6:'Fetching rpairs',
               7:'Fetching reactions - compounds links',
               8:'Fetching details on KEGG entries',
               9:'Crafting results'}
    
    _substatuses = [3,4,5,6,7,8]
    
    def __init__(self, threads=40, avoid=[], keeptrying=False,
                        queue=Queue.Queue()):
        BaseMapper.__init__(self, threads=threads, avoid=avoid,
                                    keeptrying=keeptrying, queue=queue)
        
    def getAllPathways(self):
        '''
        Get all the available pathway IDs
        '''
        kegg = KeggAPI()
        
        kegg.getIDListFromDB('pathway')
        
        for p in kegg.result:
            self.pathdet[p] = None
    
    def run(self):
        self.updateStatus()
        try:
            self.checkConnection()
        except Exception as e:
            self.sendFailure(str(e))
            return
    
        # Get pathways
        self.updateStatus()
        try:
            self.getAllPathways()
        except Exception as e:
            self.sendFailure(str(e))
            return
        
        if self.killed:
            return
        
        # Compounds
        self._maxsubstatus = len(self.pathdet)
        self.updateStatus()
        try:
            self.getPathCompounds()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Reactions for each compound
        self._maxsubstatus = len(self.compdet)
        self.updateStatus()
        try:
            self.getCompoundReacts()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Reactions
        self._maxsubstatus = len(self.pathdet)
        self.updateStatus()
        try:
            self.getPathReactions()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Related rpairs
        self._maxsubstatus = len(self.reactdet)
        self.updateStatus()
        logger.warning('Using RCLASS attribute of KEGG reactions, as the RPAIR database has now been discontinued')
        try:
            self.getReactRPairs()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Compounds for each reaction
        self._maxsubstatus = len(self.reactdet)
        self.updateStatus()
        try:
            self.getReactCompounds()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Details
        # Pathway details
        self._maxsubstatus = len(self.pathdet)
        self.updateStatus()
        try:
            self.getPathDetails()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Pathway HTML maps (!!!)
        self._maxsubstatus = len(self.pathdet)
        try:
            self.getMapsDetails()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Reaction details
        self._maxsubstatus = len(self.reactdet)
        try:
            self.getReactDetails()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Compound details
        self._maxsubstatus = len(self.compdet)
        try:
            self.getCompDetails()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # RPair details
        self._maxsubstatus = len(self.rpairdet)
        logger.warning('Using RCLASS attribute of KEGG reactions, as the RPAIR database has now been discontinued')
        try:
            self.getRPairDetails()
        except Exception as e:
            self.sendFailure(str(e))
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Prepare the output object
        self.updateStatus()
        self.result = KeggDetails()
        self.result.setDetails(None, self.reactdet,
                               self.compdet, self.pathdet, self.rpairdet)
        self.result.setLinks(pathreact=self.pathreact, 
                             pathcomp=self.pathcomp, reactcomp=self.reactcomp,
                             compreact=self.compreact,
                             rpairreact=self.rpairreact,
                             reactrpair=self.reactrpair)
        self.result.setMaps(self.pathmap)
