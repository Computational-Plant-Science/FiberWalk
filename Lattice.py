'''
FIBER WALK DEMO

This class generates the growing lattice. 
It expands the lattice around every reached point.

The demo uses:
- the networkX package (http://networkx.lanl.gov/)
- the numpy package (http://sourceforge.net/projects/numpy/)
- the scipy package (http://www.scipy.org/SciPy)
- the classes Lattice and FiberWalk written by Alexander Bucksch

Note: The graphics shown in the paper where generated with mayavi2
(http://docs.enthought.com/mayavi/mayavi/index.html)

The code is free for non-commercial use.
Please contact the author for commercial use.

Please cite the Fiber Walk Paper if you use the code for your scientific project.

-------------------------------------------------------------------------------------------
Author: Alexander Bucksch
School of Biology and Interactive computing
Georgia Institute of Technology

Mail: bucksch@gatech.edu
Web: http://www.bucksch.nl

-------------------------------------------------------------------------------------------
Copyright (c) 2012 Alexander Bucksch
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above
    copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided
    with the distribution.

  * Neither the name of the Fiber Walk Demo Developers nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

#!/usr/bin/python
import numpy as np
import networkx as nx


class Lattice:
    
    def __init__(self, dim, size=0, resetable=False):
        self.__test=5
        self.__junkSize = 500
        self.__lostNodes = []
        self.__dimension = dim
        self.__size = size
        self.__GD = nx.MultiDiGraph()
        pos = []
        for i in range(dim):
            pos.append(float(0.0))
        self.__GD.add_node(tuple(pos), consumed=False, counter=0, SA = False,Contr=1, pos=pos)
        self.expandLattice(pos)
        #self.generateLattice()
        
        
    def get_dimension(self):
        return self.__dimension

    def add_lostNode(self, node):
        self.__lostNodes.append(node)
    def get_lostNodes(self):
        return self.__lostNodes
    def is_lostNode(self,node):
        for i in self.__lostNodes:
            if i == node: return True
        return False
    
    def get_size(self):
        return self.__size


    def get_Lattice(self):
        return self.__GD

    def reset_Lattice(self):
        self.__GD.clear()
        self.__lostNodes =[]
        pos = []
        for i in range(self.dim):
            pos.append(0)
        self.expandLattice(pos)

    def set_dimension(self, value):
        self.__dim = value


    def set_size(self, value):
        self.__size = value


    def set_Lattice(self, value):
        self.__GD = value

    def del_dimension(self):
        del self.__dim


    def del_size(self):
        del self.__size


    def del_Lattice(self):
        del self.__GD
        
    # overload []
    def __getitem__(self, index):
        print "__getitem__"
        try:
            idx=list(index)
            if len(index)>self.__dimension: raise IndexError('Dimension is unequal dimensions given by the index')
            for d in range(self.__dimension):
                idx[d] = index[d]/self.__junkSize
            path = './tmp'
            for d in range(self.__dimension):
                path=path+'_'+str(idx[d])
            print self.__GD.node[0,0]['coord']   
            self.__GD = nx.read_gpickle(path)
            print self.__GD.node[0,0]['coord']
            return self.__GD.node[0,0]['coord']
        except IndexError:
            if (len(index)<self.__dimension): print 'dimension is smaller then index dimension'
            elif (len(index)>self.__dimension): print 'dimension is bigger then index dimension'
            
        #self.data[key] = item

    # overload set []
    def __setitem__(self, key, item):
        print "__setitem__"
        try:
            if len(key)>self.__dimension: raise IndexError('Dimension is unequal dimensions given by the index')
            idx=list(key)
            if len(key)>self.__dimension: raise IndexError('Dimension is unequal dimensions given by the index')
            for d in range(self.__dimension):
                idx[d] = key[d]/self.__junkSize
            path = './tmp'
            for d in range(self.__dimension):
                path=path+'_'+str(idx[d])
            print self.__GD.node[0,0]['coord']   
            self.__GD = nx.read_gpickle(path)
            print self.__GD.node[0,0]['coord']
            print 'set item is not implemented'
            #self.__GD.node[tupel(idx)][item]=value
        except IndexError:
            if (len(key)<self.__dimension): print 'dimension is smaller then index dimension'
            elif (len(key)>self.__dimension): print 'dimension is bigger then index dimension'

    
    def generateLattice(self):
        edgeCounter = 0
        maxEdges = self.__dimension * 2 * (self.__size ** self.__dimension)
        
        # A grid graph forming the lattices
        dimGrid = []
        for i in range(self.__dimension):
            dimGrid.append(self.__size)
        G = nx.grid_graph(dimGrid)
        print 'Gridgraph constructed! size = ' + str(self.__size) 
        
        # use a classic Bi-directional graph
        print 'start adding a Stairway to Heaven '
        oldcounter = 0
        for e in G.edges_iter():
            edgeCounter+=1
            ec = edgeCounter/1000
            if ec > oldcounter:
                print '.';
                oldcounter += ec
            self.__GD.add_node(e[0], consumed=False, counter=0)
            self.__GD.add_node(e[1], consumed=False, counter=0)
            if self.__GD.has_edge(e[0], e[1]) == False:
                self.__GD.add_edge(e[0], e[1], vdir=np.array(e[0]) - np.array(e[1]),SA= False,New=True)
                
            if self.__GD.has_edge(e[1], e[0]) == False:
                self.__GD.add_edge(e[1], e[0], vdir=np.array(e[1]) - np.array(e[0]),SA=False,New=True)
        print 'clearing original graph!'
        G.clear()    
    
    def expandlatticeToNeighbours(self,pos):
            for i in range(self.dim):
                posNew=np.array(pos)
                posNew2=np.array(pos)
                posNew[i] +=1.0

                vdir=np.array(posNew)-np.array(pos)
                crit =True

                if self.is_lostNode(tuple(posNew)) == False: 
                    if  self.__GD.has_node(tuple(posNew)) == False: 
                        self.__GD.add_node(tuple(posNew), consumed=False, counter=0, SA = False,Contr=1, pos=np.array(posNew))
                        
                    if  self.__GD.has_edge(tuple(pos), tuple(posNew)) == False: 
                        self.__GD.add_edge(tuple(pos),tuple(posNew),vdir=np.array(posNew)-np.array(pos), SA = False,Contr=1)
                        
                posNew2[i] -=1.0
                vdir=np.array(posNew2)-np.array(pos)
                crit =True
   
                if self.is_lostNode(tuple(posNew2)) == False: 
                    if  self.__GD.has_node(tuple(posNew2)) == False: 
                        self.__GD.add_node(tuple(posNew2), consumed=False, counter=0, SA = False,Contr=1, pos=np.array(posNew2))
               
                    if  self.__GD.has_edge(tuple(pos), tuple(posNew2)) == False: 
                        self.__GD.add_edge(tuple(pos),tuple(posNew2),vdir=np.array(posNew2)-np.array(pos), SA = False,Contr=1)
    
    def expandLatticeNoNeighbour(self,pos):
            for i in range(self.dim):
                posNew=np.array(pos)
                posNew2=np.array(pos)
                posNew[i] +=1.0
                vdir=np.array(posNew)-np.array(pos)
                crit =True
                
                if self.is_lostNode(tuple(posNew)) == False: 
                    if  self.__GD.has_node(tuple(posNew)) == False: 
                        self.__GD.add_node(tuple(posNew), consumed=False, counter=0, SA = False,Contr=1, pos=np.array(posNew))
                        
                    if  self.__GD.has_edge(tuple(pos), tuple(posNew)) == False: 
                        self.__GD.add_edge(tuple(pos),tuple(posNew),vdir=np.array(posNew)-np.array(pos), SA = False,Contr=1)
                    
                posNew2[i] -=1.0
                vdir=np.array(pos)-np.array(posNew2)
                crit =True
   
                if self.is_lostNode(tuple(posNew2)) == False: 
                        if  self.__GD.has_node(tuple(posNew2)) == False: 
                            self.__GD.add_node(tuple(posNew2), consumed=False, counter=0, SA = False,Contr=1, pos=np.array(posNew2))

                        if  self.__GD.has_edge(tuple(pos), tuple(posNew2)) == False: 
                            self.__GD.add_edge(tuple(pos),tuple(posNew2),vdir=np.array(posNew2)-np.array(pos), SA = False,Contr=1)
                        
    def expandLattice(self,pos):
        for i in range(self.dim):
            posNew=np.array(pos)
            posNew2=np.array(pos)
            posNew[i] +=1.0

            vdir=np.array(posNew)-np.array(pos)
            crit =True
           
            if self.is_lostNode(tuple(posNew)) == False: 
                if  self.__GD.has_node(tuple(posNew)) == False: 
                    self.__GD.add_node(tuple(posNew), consumed=False, counter=0, SA = False,Contr=1,pos=np.array(posNew))

                if  self.__GD.has_edge(tuple(pos), tuple(posNew)) == False: 
                    self.__GD.add_edge(tuple(pos),tuple(posNew),vdir=np.array(posNew)-np.array(pos), SA = False,Contr=1)
                
            posNew2[i] -=1.0
            vdir=np.array(pos)-np.array(posNew2)
            crit =True

            if self.is_lostNode(tuple(posNew2)) == False: 
                    if  self.__GD.has_node(tuple(posNew2)) == False: 
                        self.__GD.add_node(tuple(posNew2), consumed=False, counter=0, SA = False,Contr=1, pos=np.array(posNew2))

                    if  self.__GD.has_edge(tuple(pos), tuple(posNew2)) == False: 
                        self.__GD.add_edge(tuple(pos),tuple(posNew2),vdir=np.array(posNew2)-np.array(pos), SA = False,Contr=1)
                    
        for i in self.__GD.neighbors(tuple(pos)):
            self.expandlatticeToNeighbours(i)
            
                
    def generateDynamicLattice(self):
        
        numberOfPickels = self.__size/self.__junkSize
        numberOfPickels = np.power(numberOfPickels,self.__dimension)
        
        dimGrid = []
        for i in range(self.__dimension):
            dimGrid.append(self.__junkSize)
        #generate IDXList
       
        G = nx.grid_graph(dimGrid)
        print 'start adding a Stairway to Heaven '
        oldcounter = 0
        edgeCounter = 0
        for e in G.edges_iter():
            edgeCounter+=1
            ec = edgeCounter/1000
            if ec > oldcounter:
                print '.';
                oldcounter += ec
            self.__GD.add_node(e[0], consumed=False, counter=0)
            self.__GD.add_node(e[1], consumed=False, counter=0)
            if  self.__GD.has_edge(e[0], e[1]) == False: self.__GD.add_edge(e[0], e[1], vdir=np.array(e[0]) - np.array(e[1]),SA=False,New=True)
            if  self.__GD.has_edge(e[1], e[0]) == False: self.__GD.add_edge(e[1], e[0], vdir=np.array(e[1]) - np.array(e[0]),SA=False,New=True)
        print 'clearing original graph from memory!'
        G.clear()         
        for i in range(numberOfPickels):
            print 'storing the graph'
            idx=np.zeros(self.__dimension)
            for d in reversed(range(self.__dimension)):
                idx[d] = i/(np.power(self.__size/self.__junkSize,d))
                if d != self.__dimension-1: 
                    idx[d]=idx[d]-idx[d+1]*(self.__size/self.__junkSize)
            
            path = './tmp'
            for d in reversed(range(self.__dimension,)):
                path=path+'_'+str(int(idx[d]))
            print path
            coords=list(idx) 
            for v in self.__GD.nodes():
                for d in range(self.__dimension):
                    print v[d]
                    coords[d]=idx[d]*self.__junkSize+v[d]
                self.__GD.node[v]['coord']=tuple(coords)
                print self.__GD.node[v]
            nx.write_gpickle(self.__GD,path)
        
    dim = property(get_dimension, set_dimension, del_dimension, "dim's docstring")
    size = property(get_size, set_size, del_size, "size's docstring")
    GD = property(get_Lattice, set_Lattice, del_Lattice, "GD's docstring")

