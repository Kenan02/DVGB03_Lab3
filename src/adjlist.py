#!/usr/bin/env python3

#
# Recommended order to get started:
# 1) AdjacencyList.{add_node,node_cardinality}
# 2) AdjacencyList._add_edge, Edge.add - note: if you run in debug mode, you
# will get all edges printed when selecting 'v: view' in the menu.
# 3) AdjacencyList.find_edge, Edge.find
#

#LABB GJORD AV KENAN SAHINOVIC & OLIVER RANER

import sys
import logging

log = logging.getLogger(__name__)

from math import inf

class AdjacencyList:
    '''
    A linked-list implementation of an adjacency list that keeps its nodes and
    edges lexicographically ordered at all times.
    '''
    def __init__(self, name=None, info=None):
        '''
        Initializes a new adjacency list.  It is considered empty if no head
        node is provided.  Optionally, a node can also have associated info.
        '''
        self._name = name # head node name
        self._info = info # head node info
        if not self.head().is_empty():
            self._tail = AdjacencyList() # empty tail
            self._edges = Edge() # empty list of edges

    def is_empty(self):
        '''
        Returns true if this adjacency list is empty.
        '''
        return self._name is None

    def head(self):
        '''
        Returns the head of this adjacency list.
        '''
        return self

    def tail(self):
        '''
        Returns the tail of this adjacency list.
        '''
        return self._tail

    def cons(self, tail):
        '''
        Returns the head of this adjacency list with a newly attached tail.
        '''
        self._tail = tail
        return self.head()

    def name(self):
        '''
        Returns the node name.
        '''
        return self._name

    def info(self):
        '''
        Returns auxilirary node info.
        '''
        return self._info

    def edges(self):
        '''
        Returns the edge head.
        '''
        return self._edges

    def set_name(self, name):
        '''
        Sets the node name to `name`.

        Returns an adjacency list head.
        '''
        self._name = name
        return self.head()

    def set_info(self, info):
        '''
        Sets the auxilirary info of this node to `info`.

        Returns an adjacency list head.
        '''
        self._info = info
        return self.head()

    def set_edges(self, edges):
        '''
        Sets the edge head of this node to `edges`.

        Returns an adjacency list head.
        '''
        self._edges = edges
        return self.head()

    ###
    # Node operations
    ###
    def add_node(self, name, info=None):
        '''
        Adds a new node named `name` in lexicographical order.  If node `name`
        is a member, its info-field is updated based on `info`.

        Returns an adjacency list head.
        '''
        
        
        if self.is_empty():
            return AdjacencyList(name, info)
        
        
        #When name is smaller than the head value
        if name < self.name():
            newHead = AdjacencyList(name, info)
            return newHead.cons(self.head())
            
        #When it is larger
        else:
            self.cons(self.tail().add_node(name, info))
            
            
        return self.head()

    def delete_node(self, name):
        '''
        Deletes the node named `name` if it is a member.

        Returns an adjacency list head.
        '''
        if self.head().is_empty():
            return self.head()
        
        if name == self.name():
            return self.tail()
        
        elif name > self.name():
            self.cons(self.tail().delete_node(name))
            
        return self.head()
       

    def find_node(self, name):
        '''
        Returns True if the node named `name` is a member.
        '''
        if self.is_empty():
            return False
        if name == self.head().name():
            return True
        return self.tail().find_node(name)

    def node_cardinality(self):
        '''
        Returns the number of nodes.
        '''
        if self.is_empty():
            return 0
        else:
            return (1 + self.tail().node_cardinality())

    ###
    # Edge operations
    ###
    def add_edge(self, src, dst, weight=1):
        '''
        Adds or updates an edge from node `src` to node `dst` with a given
        weight `weight`.  If either of the two nodes are non-members, the same
        adjacency list is returned without any modification.

        Returns an adjacency list head.
        '''
        if not self.find_node(dst):
            return self.head()
        return self._add_edge(src, dst, weight)

    def _add_edge(self, src, dst, weight):
        '''
        Adds a new (or updates an existing) edge from node `src` to node `dst`,
        setting the weight to `weight`.

        Returns an adjacency list head.

        Pre: `dst` is a member of this adjacency list.
        '''
        
        if self.head().is_empty():
            return self.head()
        if src == self.name():
            self.set_edges(self.edges().add(dst, weight))
            
        elif src > self.name():
            self.tail()._add_edge(src, dst, weight)
            
            
        return self.head()
        
       

    def delete_edge(self, src, dst):
        '''
        Deletes an edge from node `src` to node `dst` if it exists.

        Returns an adjacency list head.
        '''
        if self.head().is_empty():
            return self.head()
        
        if self.name() == src:
            self.set_edges(self.edges().delete(dst))
            
        elif src > self.name():
            self.tail().delete_edge(src, dst)
            
        return self.head()

    def delete_edges(self, name):
        '''
        Deletes all edges towards the node named `name`.

        Returns an adjacency list head.
        '''
        if self.head().is_empty():
            return self.head()
        else:
            self.set_edges(self.edges().delete(name))
            self.tail().delete_edges(name)
        return self.head()

    def find_edge(self, src, dst):
        '''
        Returns True if there's an edge from node `src` to node `dst`.
        '''
        if self.head().is_empty():
            return False
        if src == self.name():
            return self.edges().find(dst)
        else:
            return self.tail().find_edge(src, dst)
        

    def edge_cardinality(self):
        '''
        Returns the number of edges.
        '''
        if self.head().is_empty():
            return 0
        else: 
            return self.edges().cardinality() + self.tail().edge_cardinality()

    def self_loops(self):
        '''
        Returns the number of loops in this adjacency list.  Note that a loop is
        defined as a node that has an edge towards itself, e.g., A->A.
        '''
        if self.head().is_empty():
            return 0
        else:
            if self.edges().find(self.name()) == True:
                return 1 + self.tail().self_loops()
            
            else:
                return self.tail().self_loops()


    def adjacency_matrix(self):
        '''
        Returns this adjacency list as an adjacency matrix.  For example,
        consider the following adjacency list where all edges have weight=1.
        
        a: a->b->c
        |
        v
        b: a->b
        |
        v
        c: c

        Then we would expect the following 3x3 adjacency matrix:

          a b c
        -+-----
        a|1 1 1
        b|1 1 *
        c|* * 1

        Where the corresponding python-matrix is:

        [ [1,1,1], [1,1,inf], [inf,inf,1] ]

        Note that inf indicates that there is no path between two nodes.  Also,
        all rows and columns are lexicographically ordered based on node names.

        Hint: depending on your solution, you may need to add a helper method
        that maps a node's name to it's numeric position in the adjacency list.
        '''
        if self.is_empty():
            return [[]]

        # In case you'd like to create an inf-initialized n x n matrix
        n = self.node_cardinality()
        matrix = [ [inf]*n for i in range(n) ]
        
        nodeList = self.getNames()
        tempNode = self.head()
        count = 0
        
        while not tempNode.is_empty():
            tempEdge = tempNode.edges()
            for i in range(n):
                if not tempEdge.is_empty():
                    if tempEdge.dst() == nodeList[i]:
                        matrix[count][i] = tempEdge.weight()
                        tempEdge = tempEdge.tail()
                        
                else:
                    break
                
            tempNode = tempNode.tail()
            count += 1
            
        return matrix
    
    
    
    #Help method
    def getNames(self):
        
        nodeList = []
        tempNode = self.head()
        
        while not tempNode.is_empty():
            nodeList.append(tempNode.name())
            tempNode = tempNode.tail()
        return nodeList
    
    

    def list_nodes(self):
        '''
        Returns a list of node names in lexicographical order.
        '''
        head, node_names = self.head(), []
        while not head.is_empty():
            node_names += [ head.name() ]
            head = head.tail()
        return node_names

    def list_edges(self):
        '''
        Returns a list of edges in lexicographical order.
        '''
        if self.head().is_empty():
            return []
        return self.head().edges().list(self.head().name()) +\
            self.tail().list_edges()
            
            
            
#HELP FUNCTIONS
    
    def ListOfNodes(self):
        '''
        pre: self is a node
        post: returns list of nodes, not the values only
        '''
        nodeHeads = []
        tmpNode = self.head()
        
        while not tmpNode.is_empty():
            nodeHeads.append(tmpNode.head())
            tmpNode = tmpNode.tail()
        return nodeHeads
    
    def ListOfEdges(self):
        '''
        pre: self is a edge
        post: returns a list of edges not just the value
        '''
        
        EdgeList = []
        tmpEdge = self.edges()
        
        while not tmpEdge.is_empty():
            EdgeList.append(tmpEdge.head())
            tmpEdge = tmpEdge.tail()
            
        return EdgeList   
    
    def getNode(self, name):
        '''
        returns node if the node with the name "name" is a member
    
        '''
        
        if self.is_empty():
            return 0
        if name == self.head().name():
            return self.head()
        return self.tail().getNode(name)
             
    
    
#HELP FUNCTIONS

class Edge:
    '''
    A linked-list implementation of edges that originate from an implicit source
    node.  Each edge has a weight and goes towards a given destination node.
    '''
    def __init__(self, dst=None, weight=1):
        '''
        Initializes a new edge sequence.  It is considered empty if no head edge
        is provided, i.e., dst is set to None.
        '''
        self._dst = dst # where is this edge's destination
        self._weight = weight # what is the weight of this edge
        if not self.head().is_empty():
            self._tail= Edge() # empty edge tail

    def is_empty(self):
        '''
        Returns true if this edge is empty.
        '''
        return self._dst is None
    
    def head(self):
        '''
        Returns the head of this edge.
        '''
        return self

    def tail(self):
        '''
        Returns the tail of this edge.
        '''
        return self._tail

    def cons(self, tail):
        '''
        Returns the head of this sequence with a newly attached tail.
        '''
        self._tail = tail
        return self.head()

    def dst(self):
        '''
        Returns the node name that this edge goes towards.
        '''
        return self._dst

    def weight(self):
        '''
        Returns the weight of this edge.
        '''
        return self._weight

    def set_dst(self, dst):
        '''
        Sets the destination of this edge to `dst`.

        Returns an edge head.
        '''
        self._dst = dst
        return self.head()

    def set_weight(self, weight):
        '''
        Sets the weight of this edge to `weight`.

        Returns an edge head.
        '''
        self._weight = weight
        return self.head()
    
    ###
    # Operations
    ###
    def add(self, dst, weight=1):
        '''
        Adds a new edge towards `dst` in lexicographical order.  If such an
        edge exists already, the associated weight-field is updated instead.

        Returns an edge head.
        '''
        
        #When the edge head is empty
        if self.head().is_empty():
            self.set_dst(dst)
            self.set_weight(weight)
            self._tail = Edge()
            return self.head()
        
        #when edge exists, only update weight
        elif self.dst() == dst:
            self.set_weight(weight)
            return self.head()
        
        #when dst is lower than edge
        elif dst < self.dst():
            new = Edge(dst, weight)
            new._tail = self.head()
            return new.head()
        
        #else if it is higher than edge
        else:
            return self.cons(self.tail().add(dst, weight))
        
       
    def delete(self, dst):
        '''
        Deletes the edge that goes towards `dst` if it exists.

        Returns an edge head.
        '''
        if self.is_empty():
            return self.head()
        if self.dst() == dst:
            return self.tail()
        return self.cons(self.tail().delete(dst))

    def find(self, dst):
        '''
        Returns True if there is an edge towards `dst` in this sequence.
        '''
        if self.dst() is None:
            return False
        elif self.dst() == dst:
            return True
        else:
            return self.tail().find(dst)

    def cardinality(self):
        '''
        Returns the number of edges in this sequence.
        '''
        if self.head().is_empty():
            return 0
        else:
            return (1 + self.tail().cardinality())
        
    def list(self, src):
        '''
        Returns a list of edges in lexicographical order, e.g., if `src`
        goes to nodes A and B, the returned list would be:
            [ (src, A), (src, B) ]
        '''
        if self.head().is_empty():
            return []
        return [(src, self.head().dst(), self.weight())] + self.tail().list(src)

if __name__ == "__main__":
    log.critical("module contains no main method")
    sys.exit(1)
