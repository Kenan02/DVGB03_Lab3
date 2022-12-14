#!/usr/bin/env python3

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
    tempMatrix = copy.deepcopy(adjlist.adjacency_matrix())
    
    #Warshall alogrithm
    for j in range(nodesize):
        for i in range(nodesize):
            if j == i:
                tempMatrix[j][i] = 0
            for k in range(nodesize):
                tempMatrix[i][k] = min(tempMatrix[i][k], (tempMatrix[i][j] + tempMatrix[j][k]))
                
                
    for i in range(nodesize):
        matrix[i] = [False if x == inf else True for x in tempMatrix[i]]
    
    
    
    return matrix

def floyd(adjlist):
    '''
    Returns an NxN matrix that contains the result of running Floyd's algorithm.

    Pre: adjlist is not empty.
    '''
    log.info("TODO: floyd()")
    return [[]]

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
    log.info("TODO: dijkstra()")
    d = [None] * adjlist.node_cardinality()
    e = [None] * adjlist.node_cardinality()
    
    q = [] #temporary queue to to sort and pick a node for loop
    s = [] #values that pop from q is stored here to order the nodes
    
    
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
    log.info("TODO: prim()")
    l = []
    c = []
    return l, c

if __name__ == "__main__":
    logging.critical("module contains no main")
    sys.exit(1)
