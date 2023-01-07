#!/usr/bin/env python3
#LABB GJORD AV KENAN SAHINOVIC & OLIVER RANER

import sys
import logging

log = logging.getLogger(__name__)

import copy #Used for deepcopy/copy
from adjlist import AdjacencyList

from math import inf

def warshall(adjlist):
    '''
    Returns an NxN matrix that contains the result of running Warshall's
    algorithm.

    Pre: adjlist is not empty.
    '''
    
    
    nodesize = adjlist.node_cardinality() #Stores number of nodes
    matrix = [[None]*nodesize]*nodesize #Matrix for result
    tempMatrix = copy.deepcopy(adjlist.adjacency_matrix()) #Temporary matrix that the algorithm uses for calculation
    
    #Warshall alogrithm
    for j in range(nodesize):
        for i in range(nodesize):
            if j == i:
                tempMatrix[j][i] = 0
            for k in range(nodesize):
                tempMatrix[i][k] = min(tempMatrix[i][k], (tempMatrix[i][j] + tempMatrix[j][k]))
                
    #Swapping from inf/int to true/false        
    for i in range(nodesize):
        matrix[i] = [False if x == inf else True for x in tempMatrix[i]]
    
    
    
    return matrix

def floyd(adjlist):
    '''
    Returns an NxN matrix that contains the result of running Floyd's algorithm.

    Pre: adjlist is not empty.
    '''
    
    nodesize = adjlist.node_cardinality()
    matrix = adjlist.adjacency_matrix() #initialize matrix with adjecency matrix
    
    tmpMatrix = copy.deepcopy(matrix) # deepcopy to make sure no values get distorted
    
    for i in range(nodesize):
        for j in range(nodesize):
            if i == j:
                tmpMatrix[i][j] = 0
            for k in range(nodesize):
                tmpMatrix[j][k] = min(tmpMatrix[j][k], (tmpMatrix[j][i] + tmpMatrix[i][k]))
                
    return tmpMatrix
    
    

def dijkstra(adjlist, start_node):
    '''
    Returns the result of running Dijkstra's algorithm as two N-length lists:
    1) distance d: here, d[i] contains the minimal cost to go from the node
    named `start_node` to the i:th node in the adjacency list.
    2) edges e: here, e[i] contains the node name that the i:th node's shortest
    path originated from.

    If the index i refers to the start node, set the associated values to None.

    Pre: start_node is a member of adjlist.

    === Example ===
    Suppose that we have the following adjacency matrix:

      a b c
    -+-----
    a|* 1 2
    b|* * 2
    c|* * *

    For start node "a", the expected output would then be:

    d: [ None, 1, 2]
    e: [ None, 'a', 'a' ]
    '''
    d = [None] * adjlist.node_cardinality()
    e = [None] * adjlist.node_cardinality()
    
    q = [] #temporary queue to sort and pick a node for loop
    s = [] #values that pop from q is stored here to keep track of which node is being evaluated or has been evaluated
    
    #Initialize
    initArr(adjlist, q, start_node)
    
    while len(q) != 0:
        q.sort(key=lambda node: node.info()[0])
        u = q.pop(0)
        
        s.append(u)
        #Get all edges from the current node then tries to relax if v has not been visited
        
        for edge in u.ListOfEdges():
            v = adjlist.getNode(edge.dst())
            if v in q:
                relax(v, u, edge.weight())
                
    #sets value in the correct spot of e and d arrays
    s.sort(key=lambda node: node.name())
    for i, node in enumerate(s):
        if node.info()[0] == 0: #Special case for the start node
            d[i] = None
            e[i] = None
        else:
            d[i] = node.info()[0]
            e[i] = node.info()[1]
            
    
    return d, e

def prim(adjlist, start_node):
    '''
    Returns the result of running Prim's algorithm as two N-length lists:
    1) lowcost l: here, l[i] contains the weight of the cheapest edge to connect
    the i:th node to the minimal spanning tree that started at `start_node`.
    2) closest c: here, c[i] contains the node name that the i:th node's
    cheapest edge orignated from. 

    If the index i refers to the start node, set the associated values to None.

    Pre: adjlist is setup as an undirected graph and start_node is a member.

    === Example ===
    Suppose that we have the following adjacency matrix:

      a b c
    -+-----
    a|* 1 3
    b|1 * 1
    c|3 1 *

    For start node "a", the expected output would then be:

    l: [ None, 1, 1]
    c: [ None, 'a', 'b' ]
    '''
    
    l = [None] * adjlist.node_cardinality() #Lowcost, weights
    c = [None] * adjlist.node_cardinality() #closest, parents
    
    q = [] #temporary queue for nodes
    
    
    #initialize arrays and values
    initArr(adjlist, q, start_node)
    
    while len(q) != 0:
        q.sort(key=lambda node: node.info()[0]) #sorts to pop the node with smallest key
        
        u = q.pop(0)
        
        for edge in (u.ListOfEdges()):
            dstNode = adjlist.getNode(edge.dst())
            if dstNode in q and edge.weight() < dstNode.info()[0]:
                dstNode.set_info([edge.weight(), u.name()]) # set new value if dst is in q and key is larger than weight
                
                
    #adds values from node.info to right place in l and c
    
    for i, node in enumerate(adjlist.ListOfNodes()):
        if node.info()[0] == 0:
            l[i] = None
            c[i] = None
        else:
            l[i] = node.info()[0]
            c[i] = node.info()[1]
            
    return l, c


#A function inspired by the help sheet on the exam on the last question
#A initialize function for algorithms that uses array of nodes and needs
#initialization of inf and 0 in array before running algorithm
def initArr(adjNode, nodelist, start):
    #pre: adjNode = root, nodelist = list of nodes, start = value of the startnode 
    #post: nodelist is an array of all the nodes in the graph, start in the array has a special value of 0 instead if inf
    #post: all nodes in the array updates info with right values
    
    
    nodelist.extend(adjNode.ListOfNodes())
    
    for node in (nodelist):
        if node.name() == start:
            node.set_info([0, None, None])
        else:
            node.set_info([inf, None, None])
            

#Function also inspired by the help sheet from the exam on the last question which compares an edge dst info 
#if (u + weight) < v then update the weight to the new and lower weight 
def relax(v, u, weight):
    
    #pre: v and u are nodes we are checking connections between and weight is an integer
    #post: updates info of v if u.info + weight is lower than the current info
    
    if v.info()[0] > u.info()[0] + weight:
        v.set_info([(u.info()[0] + weight), u.name()])
        

if __name__ == "__main__":
    logging.critical("module contains no main")
    sys.exit(1)
