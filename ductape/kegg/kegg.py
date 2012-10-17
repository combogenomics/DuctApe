#!/usr/bin/env python
"""
Kegg

Common Library

KeggAPI handles the connection to KEGG API through wsdl
KoMapper handles reactions, pathways and compounds retrieval on Ko IDs 
"""
import urllib
from ductape.common.commonthread import CommonThread
from ductape.common.utils import get_span
from ductape.kegg.web import kheader
import Queue
import logging
import os
import shutil
import threading
import time

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('ductape.kegg')

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
    def __init__(self):
        self.baseurl = 'http://www.kegg.jp/'
        self._apiurl = 'http://rest.kegg.jp/'
        self._maplink = 'http://www.kegg.jp/kegg-bin/show_pathway?'
        self.clean()
    
    def clean(self):
        self.input = None
        self.result = None

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
                res += line.rstrip().lstrip()
            elif b and line[0] != ' ':
                b = False
                return res
        
        if res == '':
            return None
        return res
    
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

    def getTitle(self, entry, retries=5):
        '''
        Get the title of a specific KEGG object
        '''
        attempts = 0
        while True:
            try:
                self.input = entry
                logger.debug('Looking for title for KEGG entry %s'%entry)
                url = self._apiurl + 'get/%s'%entry
                data = urllib.urlopen(url).read()
                self.result = [self.getEntryTag(data, 'NAME'),
                               self.getEntryTag(data, 'DEFINITION')]
                return
            except Exception, e:
                attempts += 1
                logger.debug('get failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep(2*attempts)
                if attempts >= retries:
                    logger.warning('get failed!')
                    raise Exception('get request failed')
    
    def getReactions(self, ko_ids, retries=5):
        '''
        Get the reaction IDs for a given KO list
        '''
        attempts = 0
        while True:
            try:
                self.input = ko_ids
                logger.debug('Looking for KEGG reactions from %d KO IDs'%len(ko_ids))
                url = self._apiurl + 'link/reaction/'
                for ko_id in ko_ids:
                    url += '%s+'%ko_id
                url = url.rstrip('+')
                
                data = urllib.urlopen(url).read()
                self.result = self.parseLinks(data)
                return
            except Exception, e:
                attempts += 1
                logger.debug('link (reaction) failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep(2*attempts)
                if attempts >= retries:
                    logger.warning('link (reaction) failed!')
                    raise Exception('link (reaction) request failed')
                
    def getPathways(self, re_ids, retries=5):
        '''
        Get the pathway IDs for a given reaction list
        '''
        attempts = 0
        while True:
            try:
                self.input = re_ids
                logger.debug('Looking for KEGG pathways from %d RE IDs'%len(re_ids))
                url = self._apiurl + 'link/pathway/'
                for re_id in re_ids:
                    url += '%s+'%re_id
                url = url.rstrip('+')
                
                data = urllib.urlopen(url).read()
                self.result = self.parseLinks(data)
                return
            except Exception, e:
                attempts += 1
                logger.debug('link (pathway) failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep(2*attempts)
                if attempts >= retries:
                    logger.warning('link (pathway) failed!')
                    raise Exception('link (pathway) request failed')
    
    def getReactionsByComp(self, co_ids, retries=5):
        '''
        Get the reactions IDs for a given compound list
        '''
        attempts = 0
        while True:
            try:
                self.input = co_ids
                logger.debug('Looking for KEGG reactions from %d CO IDs'%len(co_ids))
                url = self._apiurl + 'link/reaction/'
                for co_id in co_ids:
                    url += '%s+'%co_id
                url = url.rstrip('+')
                
                data = urllib.urlopen(url).read()
                self.result = self.parseLinks(data)
                return
            except Exception, e:
                attempts += 1
                logger.debug('link (reaction) failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep(2*attempts)
                if attempts >= retries:
                    logger.warning('link (reaction) failed!')
                    raise Exception('link (reaction) request failed')
    
    def getReactionsFromPath(self, path_ids, retries=5):
        '''
        Get the reaction IDs for a given pathway list
        '''
        attempts = 0
        while True:
            try:
                self.input = path_ids
                logger.debug('Looking for KEGG reactions from %d PATH IDs'%len(path_ids))
                url = self._apiurl + 'link/reaction/'
                for path_id in path_ids:
                    url += '%s+'%path_id
                url = url.rstrip('+')
                
                data = urllib.urlopen(url).read()
                self.result = self.parseLinks(data)
                return
            except Exception, e:
                attempts += 1
                logger.debug('link (reaction) failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep(2*attempts)
                if attempts >= retries:
                    logger.warning('link (reaction) failed!')
                    raise Exception('link (reaction) request failed')
    
    def getCompoundsFromReaction(self, re_ids, retries=5):
        '''
        Get the compound IDs for a given reaction list
        '''
        attempts = 0
        while True:
            try:
                self.input = re_ids
                logger.debug('Looking for KEGG compounds from %d RE IDs'%len(re_ids))
                url = self._apiurl + 'link/compound/'
                for re_id in re_ids:
                    url += '%s+'%re_id
                url = url.rstrip('+')
                
                data = urllib.urlopen(url).read()
                self.result = self.parseLinks(data)
                return
            except Exception, e:
                attempts += 1
                logger.debug('link (compound) failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep(2*attempts)
                if attempts >= retries:
                    logger.warning('link (compound) failed!')
                    raise Exception('link (compound) request failed')
    
    def getCompoundsFromPath(self, path_ids, retries=5):
        '''
        Get the compound IDs for a given pathway list
        '''
        attempts = 0
        while True:
            try:
                self.input = path_ids
                logger.debug('Looking for KEGG compounds from %d PATH IDs'%len(path_ids))
                url = self._apiurl + 'link/compound/'
                for path_id in path_ids:
                    url += '%s+'%path_id
                url = url.rstrip('+')
                
                data = urllib.urlopen(url).read()
                self.result = self.parseLinks(data)
                return
            except Exception, e:
                attempts += 1
                logger.debug('link (compound) failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep(2*attempts)
                if attempts >= retries:
                    logger.warning('link (compound) failed!')
                    raise Exception('link (compound) request failed')
    
    def getHTMLColoredPathway(self, path_id, obj_list, color_list, retries=5):
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
                url = self._maplink + path_id.lstrip('path:') + '/default%3white/'
                for i in range(len(obj_list)):
                    url += obj_list[i] + '%09' + color_list[i] + '/'
                
                sock=urllib.urlopen(url)
                self.result = sock.read()
                sock.close()
                return
            except Exception, e:
                attempts += 1
                logger.debug('show_pathway failed! Attempt %d'
                              %attempts)
                logger.debug('%s'%str(e))
                time.sleep(2*attempts)
                if attempts >= retries:
                    logger.warning('show_pathway failed!')
                    raise Exception('show_pathway request failed')

class KeggColor(object):
    '''
    Class KeggColor
    Holds the color information to be passed to MapsFetcher
    One object for each pathway
    '''
    def __init__(self, path, htmlmap= '', reactions={}, compounds={}):
        self.path = path 
        self.htmlmap = htmlmap
        self.reactions = reactions
        self.compounds = compounds
    
    def setMap(self, htmlmap):
        self.htmlmap = htmlmap
    
    def setReactions(self, reactions):
        self.reactions = reactions
        
    def setCompounds(self, compounds):
        self.compounds = compounds
        
    def getAll(self):
        '''
        Returns a tuple --> objects, color
        '''
        objs = [x for x in self.reactions]
        colors = [self.reactions[x] for x in self.reactions]
        objs += [x for x in self.compounds]
        colors += [self.compounds[x] for x in self.compounds]
        
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
        # Links
        self.koreact = None
        self.pathreact = None
        self.pathcomp = None
        self.compreact = None
        self.reactcomp = None
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
    
    def setDetails(self, ko=None, react=None, comp=None, path=None):
        self.ko = self._purgeDetails(ko)
        self.react = self._purgeDetails(react)
        self.comp = self._purgeDetails(comp)
        self.path = self._purgeDetails(path)
        
    def setLinks(self, koreact=None, pathreact=None, pathcomp=None,
                 compreact=None, reactcomp=None):
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
    def __init__(self, threads=5, queue=Queue.Queue()):
        CommonThread.__init__(self,queue)
        
        # Kegg connection
        self.handlers = []
        self._hindex = 0
        self.numThreads = threads
        
        for i in range(self.numThreads):
            obj = KeggAPI()
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
            
class BaseMapper(BaseKegg):
    def __init__(self, threads=5, avoid=[], queue=Queue.Queue()):
        BaseKegg.__init__(self, threads=threads, queue=queue)

        # Skip these
        self.avoid = avoid
        
        # Results
        self.reactdet = {}
        self.pathdet = {}
        self.pathreact = {}
        self.pathcomp = {}
        self.pathmap = {}
        self.compdet = {}
        self.reactpath = {}
        self.reactcomp = {}
        self.compreact = {}
        
        # Output
        self.result = None
        
    def getReactDetails(self):
        for piece in get_span(self.reactdet.keys(), self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += self.numThreads
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for react in piece:
                if react in self.avoid:
                    continue
                
                obj = threading.Thread(
                            target = self.handlers[piece.index(react)].getTitle,
                            args = (react,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if not handler.result:
                    continue
                self.reactdet[handler.input] = handler.result
                
    def getPathDetails(self):
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
                            target = self.handlers[piece.index(path)].getTitle,
                            args = (path,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if not handler.result:
                    continue
                self.pathdet[handler.input] = handler.result
    
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
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if not handler.result:
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
                
                obj = threading.Thread(
                                target = self.getHandler().getReactionsFromPath,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if not handler.result:
                    return
                
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
                
                obj = threading.Thread(
                                target = self.getHandler().getCompoundsFromPath,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if not handler.result:
                    return
                
                for path, comps in handler.result.iteritems():
                    if path not in self.pathcomp:
                        self.pathcomp[path] = comps
                comps = set([v for vs in handler.result.itervalues() for v in vs])
                for comp in comps:
                    if comp not in self.compdet:
                        self.compdet[comp] = None
                        
    def getCompDetails(self):
        for piece in get_span(self.compdet.keys(), self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += self.numThreads
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for comp in piece:
                obj = threading.Thread(
                            target = self.handlers[piece.index(comp)].getTitle,
                            args = (comp,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if not handler.result:
                    continue
                self.compdet[handler.input] = handler.result
    
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
                
                obj = threading.Thread(
                                target = self.getHandler().getPathways,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if not handler.result:
                    return
                
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
                
                obj = threading.Thread(
                                target = self.getHandler().getCompoundsFromReaction,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if not handler.result:
                    return
                
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
                
                obj = threading.Thread(
                                target = self.getHandler().getReactionsByComp,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if not handler.result:
                    return
                
                for comp, reacts in handler.result.iteritems():
                    if comp not in self.compreact:
                        self.compreact[comp] = reacts
                reacts = set([v for vs in handler.result.itervalues() for v in vs])
                for react in reacts:
                    if react not in self.reactdet:
                        self.reactdet[react] = None

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
               1:'Fetching reactions',
               2:'Fetching pathways',
               3:'Fetching pathways content',
               4:'Fetching reactions - compounds links',
               5:'Fetching details on KEGG entries',
               6:'Crafting results'}
    
    _substatuses = [1,2,3,4,5]
    
    def __init__(self, ko_list, threads=50, avoid=[], queue=Queue.Queue()):
        BaseMapper.__init__(self, threads=threads, avoid=avoid, queue=queue)
        # Kegg
        self.ko = ko_list
        
        # Results
        self.kodet = {}
        self.koreact = {}
    
    def getKOdet(self):
        for piece in get_span(self.ko, self.numThreads):
            if self.killed:
                logger.debug('Exiting for a kill signal')
                return
            
            self.cleanHandlers()
            self._substatus += self.numThreads
            if self._substatus > self._maxsubstatus:
                self._substatus = self._maxsubstatus
            self.updateStatus(sub=True)
            
            threads = []
            for ko in piece:
                if ko in self.avoid:
                    continue
                
                obj = threading.Thread(
                            target = self.handlers[piece.index(ko)].getTitle,
                            args = (ko,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if not handler.result:
                    continue
                self.kodet[handler.input] = handler.result
                
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
                
                obj = threading.Thread(
                                target = self.getHandler().getReactions,
                                args = (ids,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if not handler.result:
                    return
                
                for ko, reacts in handler.result.iteritems():
                    if ko not in self.koreact:
                        self.koreact[ko] = reacts
                reacts = set([v for vs in handler.result.itervalues() for v in vs])
                for react in reacts:
                    if react not in self.reactdet:
                        self.reactdet[react] = None
    
    def run(self):
        # Reactions
        self._maxsubstatus = len(self.ko)
        self.updateStatus()
        try:
            self.getReactions()
        except Exception, e:
            self.sendFailure(e.message)
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
        except Exception, e:
            self.sendFailure(e.message)
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
        except Exception, e:
            self.sendFailure(e.message)
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # 2. Compounds
        self._maxsubstatus = len(self.pathdet)
        try:
            self.getPathCompounds()
        except Exception, e:
            self.sendFailure(e.message)
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
        except Exception, e:
            self.sendFailure(e.message)
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Reactions for each compound
        self._maxsubstatus = len(self.compdet)
        try:
            self.getCompoundReacts()
        except Exception, e:
            self.sendFailure(e.message)
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
        except Exception, e:
            self.sendFailure(e.message)
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Pathway details
        self._maxsubstatus = len(self.pathdet)
        try:
            self.getPathDetails()
        except Exception, e:
            self.sendFailure(e.message)
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Pathway HTML maps (!!!)
        self._maxsubstatus = len(self.pathdet)
        try:
            self.getMapsDetails()
        except Exception, e:
            self.sendFailure(e.message)
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Reaction details
        self._maxsubstatus = len(self.reactdet)
        try:
            self.getReactDetails()
        except Exception, e:
            self.sendFailure(e.message)
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Compound details
        self._maxsubstatus = len(self.compdet)
        try:
            self.getCompDetails()
        except Exception, e:
            self.sendFailure(e.message)
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Prepare the output object
        self.updateStatus()
        self.result = KeggDetails()
        self.result.setDetails(self.kodet, self.reactdet,
                               self.compdet, self.pathdet)
        self.result.setLinks(koreact=self.koreact, pathreact=self.pathreact, 
                             pathcomp=self.pathcomp, reactcomp=self.reactcomp,
                             compreact=self.compreact)
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
               1:'Fetching reactions',
               2:'Fetching pathways',
               3:'Fetching pathways content',
               4:'Fetching reactions - compounds links',
               5:'Fetching details on KEGG entries',
               6:'Crafting results'}
    
    _substatuses = [1,2,3,4,5]
    
    def __init__(self, co_list, threads=50, avoid=[], queue=Queue.Queue()):
        BaseMapper.__init__(self, threads=threads, avoid=avoid, queue=queue)
        # Kegg
        self.co = co_list
        
        # Results
        self.comppath = {}
    
    def run(self):
        # Reactions
        for co_id in self.co:
            self.compdet[co_id] = None
        self._maxsubstatus = len(self.compdet)
        self.updateStatus()
        try:
            self.getCompoundReacts()
        except Exception, e:
            self.sendFailure(e.message)
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
        except Exception, e:
            self.sendFailure(e.message)
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
        except Exception, e:
            self.sendFailure(e.message)
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # 2. Compounds
        self._maxsubstatus = len(self.pathdet)
        try:
            self.getPathCompounds()
        except Exception, e:
            self.sendFailure(e.message)
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
        except Exception, e:
            self.sendFailure(e.message)
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Reactions for each compund
        self._maxsubstatus = len(self.compdet)
        try:
            self.getCompoundReacts()
        except Exception, e:
            self.sendFailure(e.message)
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
        except Exception, e:
            self.sendFailure(e.message)
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Pathway HTML maps (!!!)
        self._maxsubstatus = len(self.pathdet)
        try:
            self.getMapsDetails()
        except Exception, e:
            self.sendFailure(e.message)
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Reaction details
        self._maxsubstatus = len(self.reactdet)
        try:
            self.getReactDetails()
        except Exception, e:
            self.sendFailure(e.message)
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Compound details
        self._maxsubstatus = len(self.compdet)
        try:
            self.getCompDetails()
        except Exception, e:
            self.sendFailure(e.message)
            return
        self.cleanHandlers()
        self.resetSubStatus()
        
        if self.killed:
            return
        
        # Prepare the output object
        self.updateStatus()
        self.result = KeggDetails()
        self.result.setDetails(react=self.reactdet,
                               comp=self.compdet, path=self.pathdet)
        self.result.setLinks(pathreact=self.pathreact, 
                             pathcomp=self.pathcomp, compreact=self.compreact,
                             reactcomp=self.reactcomp)
        self.result.setMaps(self.pathmap)

class MapsFetcher(BaseKegg):
    '''
    Class MapsFetcher
    Download colored Kegg maps (png or URLs)
    Input: color_objs (KeggColor list), picture, htmls, urls, prefix, legend (file)
    Output: tuple(list of png filenames, list of HTML files, list of URLs)
    '''
    
    _statusDesc = {0:'Not started',
               1:'Making room',
               2:'Fetching maps (pictures)',
               3:'Generating interactive web pages'}
    
    _substatuses = [2]
    
    def __init__(self, color_objs, pictures=True, html=True, prefix='', 
                 legend=None, threads=20, queue=Queue.Queue()):
        BaseKegg.__init__(self, threads=threads, queue=queue)
        
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
                objs,colors = kmap.getAll()
                
                obj = threading.Thread(
                        target = self.handlers[piece.index(kmap)].getHTMLColoredPathway,
                        args = (path,objs,colors,))
                obj.start()
                threads.append(obj)
            time.sleep(0.01)
            while len(threads) > 0:
                for thread in threads:
                    if not thread.isAlive():
                        threads.remove(thread)
            for handler in self.handlers:
                if not handler.result:
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

                        sock=urllib.urlopen(urlimage)
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
                    line = './path:'.join([s[0]] + [line1])
                
                newhtml.append(line)
                
            fOut.write('%s\n'%'\n'.join(newhtml))
            fOut.write('<div id="poplay" class="poplay" />\n</body>\n</html>')
            fOut.close()
            
            self.webpages.append(fname)
    
    def run(self):
        self.updateStatus()
        self.makeRoom()
        
        if self.killed:
            return
        
        if self.pictures:
            self._maxsubstatus = len(self.colors)
            self.updateStatus()
            try:
                self.getMaps()
            except Exception, e:
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
            except Exception, e:
                self.sendFailure(e)
                return
        else:
            self.updateStatus(send=False)

class KeggNet(CommonThread):
    pass