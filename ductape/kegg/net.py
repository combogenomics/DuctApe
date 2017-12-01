#!/usr/bin/env python
"""
Net

Kegg Library

Networks made using Kegg data
"""
import logging
import networkx as nx
import numpy as np
# Nodes color handling
import matplotlib.colors as pltcls
import matplotlib.pyplot as plt

__author__ = "Marco Galardini"

################################################################################
# Log setup

logger = logging.getLogger('ductape.net')

################################################################################
# Classes

class Compound(object):
    '''
    Kegg compound, enconded as a node
    '''
    def __init__(self, co_id, name, weight=None, vmax=9):
        self.co_id = co_id
        self.name = name
        
        self.vmax = vmax

        if weight is not None:
            self.weight = weight
            
    def getColor(self):
        '''
        Transform the node weight into an hex color
        '''
        phenomeNorm = pltcls.Normalize(vmin=0, vmax=self.vmax)
        if self.weight:
            return pltcls.rgb2hex( 
                                plt.get_cmap('RdYlGn')(phenomeNorm(self.weight))
                                )
        
        return '#808080'
                    
class Reaction(object):
    '''
    Kegg reaction, enconded as an edge
    '''
    def __init__(self, re_id, co1, co2, name, weight=None):
        self.re_id = re_id
        self.co1 = co1
        self.co2 = co2
        self.name = name
        
        if weight:
            self.weight = weight

class MetabolicNet(object):
    '''
    A metabolic network as a networkX graph
    
    Nodes: kegg compounds (w co_id, name; w or w/o weight)
    Edges: kegg reactions (w co1, co2, re_id, name; w or w/o weight)
    
    nodes weight indicate the activity index
    edges weight indicate the copy number 
    '''
    def __init__(self, nodes=None, edges=None, name='MetNet'):
        self.name = name
        
        self.net = nx.Graph()
        
        if nodes:
            for n in nodes:
                self.net.add_node(n.co_id, name=n.name)
                if hasattr(n, 'weight'):
                    self.net.node[n.co_id]['weight'] = n.weight
                    self.net.node[n.co_id]['graphics'] = {'fill': n.getColor()}
        
        if edges:
            for e in edges:
                self.net.add_edge(e.co1, e.co2, reid=e.re_id, name=e.name)
                if hasattr(e, 'weight'):
                    self.net.adj[e.co1][e.co2]['weight'] = e.weight
                    
    def hasNodesWeight(self):
        '''
        At least one node has weight?
        '''
        for co in self.net.nodes():
            try:
                self.net.node[co]['weight']
                return True
            except:pass
            
        return False
    
    def hasEdgesWeight(self):
        '''
        At least one edge has weight?
        '''
        for e in self.net.edges():
            if 'weight' in self.net.adj[e[0]][e[1]]:
                return True

        return False
    
    def removeSingletons(self):
        '''
        Remove nodes with degree 0
        '''
        to_remove = filter(lambda x: self.net.degree()[x] == 0,
                           self.net.nodes())
        self.net.remove_nodes_from(to_remove)
            
    def setNet(self, net):
        '''
        Use an external networkx graph
        '''
        self.net = net
                    
    def addNodes(self, nodes):
        '''
        Takes a compounds iterable and adds (or updates) the nodes
        w co_id, name; w or w/o weight attributes
        nodes weight indicate the activity index
        '''
        for n in nodes:
            self.net.add_node(n.co_id, name=n.name)
            if hasattr(n, 'weight'):
                self.net.node[n.co_id]['weight'] = n.weight
                self.net.node[n.co_id]['graphics'] = {'fill': n.getColor()}
                
    def __len__(self):
        '''
        Returns the number of reactions
        '''
        i = 0
        for e in self.net.edges():
            if 'weight' in self.net.adj[e[0]][e[1]]:
                i += self.net.adj[e[0]][e[1]]['weight']
            else:
                i += 1
                
        return i
    
    def getDistinctReactions(self):
        '''
        Returns the distinct reaction IDs of this network
        '''
        re_id = set()
        for e in self.net.edges():
            re_id.add(self.net.adj[e[0]][e[1]]['reid'])
        
        return re_id
    
    def mean(self):
        '''
        Get the mean compounds activity (nodes weight)  
        '''
        weights = []
        
        for co in self.net.nodes():
            try:
                weights.append(self.net.node[co]['weight'])
            except:pass
            
        if len(weights) == 0:
            return np.nan
        
        return np.array(weights).mean()
    
    def std(self):
        '''
        Get the stddev compounds activity (nodes weight)  
        '''
        weights = []
        
        for co in self.net.nodes():
            try:
                weights.append(self.net.node[co]['weight'])
            except:pass
            
        if len(weights) == 0:
            return np.nan
        
        return np.array(weights).std()

    def getComponents(self):
        return len(list(nx.connected_components(self.net)))
    
    def getComponentsSizes(self):
        return [len(x) for x in nx.connected_components(self.net)]
    
    def getComponentsMean(self):
        return np.array([len(x)
                         for x in nx.connected_components(self.net)]).mean()
                         
    def getComponentsStd(self):
        return np.array([len(x)
                         for x in nx.connected_components(self.net)]).std()
