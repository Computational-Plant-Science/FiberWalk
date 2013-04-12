'''
FIBER WALK DEMO

The Fiber Walk class to run a Fiber Walk simulation.

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

import Lattice as L
import os
import pickle
import numpy as np
import networkx as nx
from scipy import polyfit, polyval, sqrt

class FiberWalk():
    
    def __init__(self, dim, numberOfObjects=1, numberOfContractions=1, minLength=0):
        #The Lattice object
        self.__lattice = L.Lattice( dim )
        #The walk itself as a networkX object
        self.__Walk = nx.Graph()
        
        # dimension of the walk
        self.__dim = dim
        self.__minLength=minLength
        #number of contractions per step
        self.__numberOfContractions = numberOfContractions
        self.__WalkPosArray = []
        self.__distancesToOrigin = []
        self.__allDistancesToOrigin = []
        self.__validNeighbors = []
        self.__contractNeighbors = []
        self.__SA = []
        self.__SAPath = []
        self.__Contr = []
        self.__pathLables = []
        
        #number of Walks to be computed
        self.__numberOfObjects = numberOfObjects
        # Distances to the origin for all walks computed
        self.__allWalksDist = []
        self.__allValidNeighbors = []
        self.__allContractNeighbors = []
        self.__allSA = []
        self.__allSAPaths = []
        self.__allContr = []
        self.__allWalksSA = []
        self.__allWalksSAPaths = []
        self.__allWalksContr = []
        
        # Norm used to calculate the distance from the seed
        self._norm = 2
        self.__saCount = 0
        self.__saPathCount = 0
        self.__ContrCount = 0
        #checking for directories to enable restart
        try:
            os.stat('./tmp/')
        except:
            os.makedirs('./tmp/')
    #get and set methods to access stats data
    def get_DeathCountArr(self):
        return self.__DeathCountArr    
    def get_lattice(self):
        return self.__lattice
    def get_allWalksContr(self):
        return self.__allWalksContr 
    def get_allWalksDist(self):
        return self.__allWalksDist
    def get_allWalksSA(self):
        return self.__allWalksSA
    def get_allWalksSAPaths(self):
        return self.__allWalksSAPaths
    def get_ContrCount(self):
        return self.__ContrCount 
    def get_SAPathCount(self):
        return self.__saPathCount   
    def get_SACount(self):
        return self.__saCount
    def get_allSA(self):
        return self.__allSA
    def get_allSAPaths(self):
        return self.__allSAPaths
    def get_allContr(self):
        return self.__allContr
    def get_size(self):
        return self.__size
    def get_dim(self):
        return self.__dim
    def get_walk(self):
        return self.__Walk
    def get_distances_to_origin(self):
        return self.__distancesToOrigin
    def get_all_distances_to_origin(self):
        return self.__allDistancesToOrigin
    def get_valid_neighbors(self):
        return self.__validNeighbors
    def get_all_valid_neighbors(self):
        return self.__allValidNeighbors
    def get_contract_neighbors(self):
        return self.__contractNeighbors
    def get_all_contract_neighbors(self):
        return self.__allContractNeighbors
    def get_path_lables(self):
        return self.__pathLables
    def set_lattice(self, value):
        self.__lattice = value
    def set_size(self, value):
        self.__size = value
    def set_dim(self, value):
        self.__dim = value
    def set_walk(self, value):
        self.__Walk = value
    def set_distances_to_origin(self, value):
        self.__distancesToOrigin = value
    def set_all_distances_to_origin(self, value):
        self.__allDistancesToOrigin = value
    def set_valid_neighbors(self, value):
        self.__validNeighbors = value
    def set_contract_neighbors(self, value):
        self.__contractNeighbors = value
    def set_path_lables(self, value):
        self.__pathLables = value
    #destructing methods
    def del_lattice(self):
        del self.__lattice
    def del_size(self):
        del self.__size
    def del_dim(self):
        del self.__dim
    def del_walk(self):
        del self.__Walk
    def del_distances_to_origin(self):
        del self.__distancesToOrigin
    def del_all_distances_to_origin(self):
        del self.__allDistancesToOrigin
    def del_valid_neighbors(self):
        del self.__validNeighbors
    def del_contract_neighbors(self):
        del self.__contractNeighbors
    def del_path_lables(self):
        del self.__pathLables
    #functions to reset arrays
    def resetArrays(self, latticeReset=False):
        if latticeReset == True:
            self.__lattice.reset_Lattice()
        self.__Walk.clear()
        self.__WalkPosArray = []
        self.__distancesToOrigin = []
        self.__validNeighbors = []
        self.__contractNeighbors = []
        self.__pathLables = []
    def resetAllArrays(self):
        self.__allDistancesToOrigin = []
        self.__allValidNeighbors = []
        self.__allContractNeighbors = []
        self.__allSA = []
        self.__allSAPaths = []
        self.__allContr = []
    
    def defineSeed(self, pos='center'):
        ret = np.zeros(self.__dim)
        if self.__lattice != 0:
            n = self.__lattice.get_Lattice().node[tuple(list(ret))]
            #mark the seed as consumed
            n['consumed'] = True
            n['counter'] = 0
            n['pos'] = ret
            n['SA'] = False
            n['Contr'] = 1
        return tuple(list(ret))    
    def selectWalkEdgeLeadingLabelSelfAvoidContractSq(self, newPos, oldPos, newWalk):
        stop = False        
        neighbors = self.__lattice.get_Lattice().neighbors(newPos)
        backEdge = self.__lattice.get_Lattice()[newPos][oldPos][0]['vdir']
        # COLLECT NEIGHBOURS AND EXCLUDE BACK EDGES AND SELFAVOIDING EDGES
        validNeighbors = []
        for i in neighbors:
            if self.__lattice.get_Lattice()[newPos][i][0]['SA'] == False:
                if self.__Walk.has_edge(newPos, i) == True:
                    print 'BUG !!!!!!'
                e = self.__lattice.get_Lattice()[newPos][i][0]['vdir']
                #check if the walk does not form a loop -> self avoidance
                if self.__lattice.get_Lattice().node[i]['consumed'] == False:
                        testBack = np.zeros_like(backEdge)
                        for k in range(len(backEdge)):
                            testBack[k] = backEdge[k] + e[k]
                            testBack[k] = np.abs(testBack[k])
                        testBack = np.sum(testBack)
                        if testBack <= self.__dim:
                            validNeighbors.append(i)
                else: 
                    if self.__lattice.get_Lattice()[newPos][i][0]['SA'] == False:
                        if self.__Walk.has_edge(newPos, i) == False:
                            self.__saPathCount +=1
                        self.__lattice.get_Lattice()[newPos][i][0]['SA'] = True
                        self.__lattice.get_Lattice()[i][newPos][0]['SA'] = True
                    
            elif newWalk == False:
                if self.__Walk.has_edge(i, newPos) == True:
                    print 'BUG !!!!!!'
                            
        if len(validNeighbors) == 0 :
            return  newPos, newPos, True, False
            
            
        oldPos = newPos
        
        rIdx = np.random.randint(0, len(validNeighbors));
        newPos = validNeighbors[rIdx]
        
        #check for valid contractions -> Just Degug code
        for i in range(self.__numberOfContractions):
                self.get_lattice().expandLattice(newPos)
                contrOK=self.contract(newPos,oldPos)
                if contrOK==False:
                    return  newPos, newPos, True, False
                
        self.get_lattice().expandLattice(newPos)
        if len(validNeighbors) != 0:
            #get the path label
            e = self.__lattice.get_Lattice()[oldPos][newPos][0]['vdir']
            #store the path lables
            self.__pathLables.append(e)
            # add edge to the walk
            if self.__Walk.has_edge(oldPos, newPos) == False: 
                self.__Walk.add_node(newPos, counter=1, pos=None,trapped=False)
                self.__Walk.add_edge(newPos, oldPos, counter=1)
                self.__lattice.get_Lattice()[newPos][oldPos][0]['SA'] = True
                self.__lattice.get_Lattice()[oldPos][newPos][0]['SA'] = True
                self.__lattice.get_Lattice()[newPos][oldPos][0]['consumed'] = True
                self.__lattice.get_Lattice()[oldPos][newPos][0]['consumed'] = True
                newWalk = True
            else:
                self.__Walk[oldPos][newPos]['counter'] += 1
                self.__Walk.node[newPos]['counter'] += 1
                newWalk = False
           
            
            self.__validNeighbors.append(len(validNeighbors))
            #Set vertex as consumed
            self.__lattice.get_Lattice().node[newPos]['consumed'] = True;
            
                
        else:  
            stop = True

        return newPos, oldPos, stop, newWalk        
    def initSeedContraction(self, newPos, oldPos):
        # Contract the Seed
        newWalk = False
        validNeighbors = []
        pathIdx = []
        contractionNeighbors = []
        #Consume the seed ->Don't consume it, otherwise it is not possible to walk out there
        self.__lattice.get_Lattice().node[oldPos]['consumed'] = True
        #contract the edges of the newly reached vertex 
        neighbors = self.__lattice.get_Lattice().neighbors(tuple(oldPos))
        for i in neighbors:
            if self.__Walk.has_edge(oldPos, i):
                pathIdx.append(i)
            elif self.__lattice.get_Lattice().node[i]['consumed'] == False:
                validNeighbors.append(i)
            elif self.__lattice.get_Lattice().node[i]['consumed'] == True:
                if self.__Walk.has_edge(oldPos, i) == False: 
                    self.__saPathCount += 1
                self.__lattice.get_Lattice()[newPos][i][0]['SA'] = True
                self.__lattice.get_Lattice()[i][newPos][0]['SA'] = True
                
        if len(validNeighbors) > 0: 
            rIdx = np.random.randint(0, len(validNeighbors));
            newPos = validNeighbors[rIdx]
       
            #make array of S to increment
            newPos = np.array(newPos)
            newPos = tuple(list(newPos)) 
            if self.__lattice.get_Lattice().node[newPos]['consumed'] == False: newWalk = True
            else: newWalk = False
            self.get_lattice().expandLattice(newPos)
            e = np.array(oldPos) - np.array(newPos)
            neighbors = self.__lattice.get_Lattice().neighbors(newPos)
            
            self.__lattice.get_Lattice().node[newPos]['consumed'] = True
    
            if self.__Walk.has_edge(oldPos, newPos) == False: 
                self.__Walk.add_node(newPos, counter=1, pos=None,trapped=False)
                self.__Walk.add_edge(oldPos, newPos, counter=1)
                self.__lattice.get_Lattice()[newPos][oldPos][0]['SA'] = True
                self.__lattice.get_Lattice()[oldPos][newPos][0]['SA'] = True
                self.__lattice.get_Lattice()[newPos][oldPos][0]['consumed'] = True
                self.__lattice.get_Lattice()[oldPos][newPos][0]['consumed'] = True
                newWalk = True
            else:
                self.__Walk[oldPos][newPos]['counter'] += 1
                self.__Walk.node[newPos]['counter'] += 1
                newWalk = False
           
            self.__pathLables.append(e)
            self.__contractNeighbors.append(len(contractionNeighbors))
            # all neighbors are valid in the init procedure
            self.__validNeighbors.append(len(neighbors))
            stop=False
        else:
            stop = True
            
            
        return newPos, oldPos, stop, newWalk
    #debug function to test if a contraction was valid
    def contractionIsValid(self, newPos, i, j):
        ret = True
        
        #trivial case
        if j == newPos:
            ret = False
        elif i == newPos:
            ret = False
        elif i == j:
            ret = False
        elif self.__lattice.get_Lattice().node[i]['consumed'] == True:
            if (self.__Walk.has_edge(newPos, i) == False):
                    self.__lattice.get_Lattice()[newPos][i][0]['SA'] = True
                    self.__lattice.get_Lattice()[i][newPos][0]['SA'] = True
                    ret = False
        #detect a possible triangle
        elif self.__lattice.get_Lattice().has_edge(newPos, j):
                if self.__lattice.get_Lattice().has_edge(newPos, i):
                    l1 = self.__lattice.get_Lattice()[newPos][i][0]['vdir']
                    l2 = self.__lattice.get_Lattice()[newPos][j][0]['vdir']
                    #check label equality
                    testEq = True
                    for k in range(len(l1)):
                        if l1[k] != l2[k]:
                            testEq = False
                            break
                    if testEq == True:
                        ret = False
                        if self.__lattice.get_Lattice()[newPos][i][0]['SA'] == False:
                            
                                self.__lattice.get_Lattice()[newPos][i][0]['SA'] = True
                                self.__lattice.get_Lattice()[i][newPos][0]['SA'] = True
                                if self.__Walk.has_edge(newPos, i) == False:self.__saCount += 1
                        if self.__lattice.get_Lattice()[newPos][j][0]['SA'] == False:
                            
                                self.__lattice.get_Lattice()[newPos][j][0]['SA'] = True
                                self.__lattice.get_Lattice()[j][newPos][0]['SA'] = True
                                if self.__Walk.has_edge(newPos, j) == False:  self.__saCount += 1
                        
                    else:
                        ret = True
                
        if ret==True:
            if self.__lattice.get_Lattice().has_edge(i, j):
                if self.__lattice.get_Lattice().has_edge(i, newPos):
                    l1 = self.__lattice.get_Lattice()[i][newPos][0]['vdir']
                    l2 = self.__lattice.get_Lattice()[i][j][0]['vdir']
                    testEq = True
                    for k in range(len(l1)):
                        if l1[k] != l2[k]:
                            testEq = False
                            break
                    if testEq == True:
                        ret = False
                        if self.__lattice.get_Lattice()[j][i][0]['SA'] == False:
                                self.__lattice.get_Lattice()[j][i][0]['SA'] = True
                                self.__lattice.get_Lattice()[i][j][0]['SA'] = True
                                self.__saCount += 1
                        if self.__lattice.get_Lattice()[i][j][0]['SA'] == False:
                                self.__lattice.get_Lattice()[i][j][0]['SA'] = True
                                self.__lattice.get_Lattice()[j][i][0]['SA'] = True
                                self.__saCount += 1
            else:
                ret = True
        if ret==True:
            if self.__lattice.get_Lattice().has_edge(j, newPos):
                l1 = self.__lattice.get_Lattice()[newPos][j][0]['vdir']
                l2 = self.__lattice.get_Lattice()[i][j][0]['vdir']
                testEq = True
                for k in range(len(l1)):
                    if l1[k] != l2[k]:
                        testEq = False
                        break
                if testEq == True:
                    ret = False
                    if self.__lattice.get_Lattice()[j][i][0]['SA'] == False:
                            self.__lattice.get_Lattice()[j][i][0]['SA'] = True
                            self.__lattice.get_Lattice()[i][j][0]['SA'] = True
                            if self.__Walk.has_edge(newPos, i) == False: self.__saCount += 1
                    if self.__lattice.get_Lattice()[j][i][0]['SA'] == False:
                            self.__lattice.get_Lattice()[j][i][0]['SA'] = True
                            self.__lattice.get_Lattice()[i][j][0]['SA'] = True
                            if self.__Walk.has_edge(newPos, i) == False: self.__saCount += 1
        else:
            ret = True 
        return ret
    
    
    # funtion to perform the contraction     
    def contract(self, newPos, oldPos):
        # Start Contraction

        contractionNeighbors = []
        # get the neighbors of the random selected vertex
        neighbors = self.__lattice.get_Lattice().neighbors(newPos)
        e= self.__lattice.get_Lattice()[oldPos][newPos][0]['vdir']
        # collect all neighbors that are possible to contract
        for i in neighbors:
            e2 = self.__lattice.get_Lattice()[newPos][i][0]['vdir']
            testEq = True
#            check the back edge of the walk
            for n in range(len(e)):
                if e[n]!=e2[n]:
                    testEq=False
#            if it is not the back edge
            if testEq==False:
                #check for self avoidance
                if self.__lattice.get_Lattice().node[i]['consumed'] == False:
                    contractionNeighbors.append(i)
                    self.__ContrCount += 1
                #check if there is space otherwise move the walk to the opposite position
                elif self.__lattice.get_Lattice()[newPos][i][0]['Contr']<self.__lattice.get_Lattice().node[newPos]['counter']:
                    return False

            
        for i in contractionNeighbors:
            n2 = self.__lattice.get_Lattice().neighbors(i)
            l1 = self.__lattice.get_Lattice()[newPos][i][0]['vdir']
            self.get_lattice().expandLattice(i)
            
            for j in n2:
                if j!=newPos:
                    l2 = self.__lattice.get_Lattice()[i][j][0]['vdir']
                    l = l1 + l2
                    for a in range(len(l)):
                        if l[a] == 0:
                            l[a] = 0
                        if l[a] < 0:
                            l[a] = -1
                        if l[a] > 0:
                            l[a] = 1
                    
                    
                    if self.contractionIsValid(newPos, i, j) == True:
                          
                        # Note: This will make every edge of the walk a non contractable edge
                        c = self.__lattice.get_Lattice()[i][j][0]['Contr']
                        if self.__lattice.get_Lattice().has_edge(newPos, j) == False:
                            self.__lattice.get_Lattice().add_edge((newPos), j, vdir=l, SA=False, Contr=c + 1)
                        if self.__lattice.get_Lattice().has_edge(j, newPos) == False:
                            self.__lattice.get_Lattice().add_edge(j, (newPos), vdir=(l * (-1)), SA=False, Contr=c + 1)
                        if self.__lattice.get_Lattice().has_edge(i, j): 
                            self.__lattice.get_Lattice().remove_edge(i, j)
                        if self.__lattice.get_Lattice().has_edge(j, i):
                            self.__lattice.get_Lattice().remove_edge(j, i)
                        self.__lattice.get_Lattice().node[newPos]['pos'] = None
                        self.__lattice.get_Lattice().node[newPos]['Contr'] += 1
                        self.__lattice.get_Lattice().node[j]['pos'] = None
            self.__lattice.get_Lattice().remove_node(i)
            self.__lattice.add_lostNode(i)
            
        return True
    
            
#The function to call for executing a Fiber Walk simulation
    def FiberWalk(self, numberOfSteps=100, frequency=2, position='center'):

        #init all Walks
        
        for n in range(self.__numberOfObjects):
            restartCount=0
            stepPossible = True
            #self.resetAllArrays()
            self.resetArrays(True)
            self.__ContrCount = 0
            self.__saPathCount = 0
            self.__saCount = 0
            seedPos = self.defineSeed(position)
            self.__Walk.add_node(seedPos,pos=(0.,0.,0.), counter=1)
            # arrays to collect statistics
            self.__WalkPosArray.append([seedPos, seedPos, False])
            self.__allDistancesToOrigin.append([0.0000001])
            self.__allValidNeighbors.append([0])
            self.__allContractNeighbors.append([0])
            self.__allSA.append([0])
            self.__allSAPaths.append([0])
            self.__allContr.append([0])   
            restartActive=False  
            stepCountTmp=len(self.__distancesToOrigin)
            stepCountMax=0
            oldStepcount=0  
            #just one step downwards to init the walk
            #some init for the run
            newWalk = self.__WalkPosArray[0][2]
            startPos = self.__WalkPosArray[0][0]
            self.__distancesToOrigin = self.__allDistancesToOrigin[n]
            self.__validNeighbors = self.__allValidNeighbors[n]
            self.__contractNeighbors = self.__allContractNeighbors[n]
            self.__SA = self.__allSA[n]
            self.__SAPath = self.__allSAPaths[n]
            self.__Contr = self.__allContr[n]  
            self.__saCount = self.__SA[len(self.__SA) - 1]
            self.__saPathCount = self.__SAPath[len(self.__SAPath) - 1]
            self.__ContrCount = self.__Contr[len(self.__Contr) - 1]
            oldPos = startPos
            
            while stepPossible == True:
               stepPossible = False
               
               #make array of currentPos to increment
               newPos = self.__WalkPosArray[0][1]

               if restartActive==False:
                   print 'step: '+str(stepCountTmp)+' of Walk Nr. '+str(n+1)+'/'+str(self.__numberOfObjects) 
                   stepCountTmp+=1
                   if stepCountMax-stepCountTmp > 20:
                       stepCountMax=stepCountTmp
               restartActive=False
               # just the first step
               if newPos == seedPos: 
                   newPos, oldPos, stop, newWalk = self.initSeedContraction(newPos, oldPos)
                   if len(self.__distancesToOrigin) < numberOfSteps and stop==False: stepPossible = True
                   else: stepPossible = False
                   if stop==True: 
                       stepPossible = False
                       restartActive=True
                   self.__WalkPosArray[0] = [oldPos, newPos, newWalk]
               # if step is not trapped, then oldPos and newPos are different
               elif oldPos != newPos:
                   #save current state to recover it in case of a restart
                   nx.write_gpickle(self.__Walk, './tmp/wtmp'+str(stepCountTmp))
                   nx.write_gpickle(self.__lattice, './tmp/ltmp'+str(stepCountTmp))
                   outfile = open('./tmp/WPAtmp'+str(stepCountTmp), 'wb')
                   pickle.dump(self.__WalkPosArray[0], outfile)
                   outfile.close()
                   outfile = open('./tmp/nptmp'+str(stepCountTmp), 'wb')
                   pickle.dump(newPos, outfile)
                   outfile.close()
                   outfile = open('./tmp/optmp'+str(stepCountTmp), 'wb')
                   pickle.dump(oldPos, outfile)
                   outfile.close()
                   outfile = open('./tmp/nwtmp'+str(stepCountTmp), 'wb')
                   pickle.dump(newWalk, outfile)
                   outfile.close()
                   outfile = open('./tmp/dtotmp'+str(stepCountTmp), 'wb')
                   pickle.dump(self.__distancesToOrigin,outfile)
                   outfile.close()
                   outfile = open('./tmp/satmp'+str(stepCountTmp), 'wb')
                   pickle.dump(self.__SA,outfile)
                   outfile.close()
                   outfile = open('./tmp/saptmp'+str(stepCountTmp), 'wb')
                   pickle.dump(self.__SAPath,outfile)
                   outfile.close()
                   outfile = open('./tmp/ctmp'+str(stepCountTmp), 'wb')
                   pickle.dump(self.__Contr,outfile)
                   outfile.close()
                   newPos, oldPos, stop, newWalk = self.selectWalkEdgeLeadingLabelSelfAvoidContractSq(newPos, oldPos, newWalk)
               elif oldPos == newPos: 
                   print 'Walk reached stopping configuration at step '+str(stepCountTmp)
                   #restart condition
                   if len(self.__distancesToOrigin) < numberOfSteps+1:
                       if stepCountTmp > 2:
                           if restartCount<numberOfSteps:
                               print 'restart the Fiber Walk'
                               if oldStepcount>=stepCountTmp: 
                                   stepCountTmp=stepCountMax-1
                                   #preventing getting negative step numbers. Step 1 is never trapped.
                                   if stepCountTmp==1: 
                                       stepCountTmp=2
                               elif oldStepcount<stepCountTmp: 
                                   oldStepcount=stepCountTmp
                                   #set the step counter one step back
                                   #Note: This occurs until the walk is longer then the first trapping
                                   stepCountTmp-=1
                                   if stepCountTmp==1: 
                                       stepCountTmp=2
                               #recover previous state
                               self.__Walk=nx.read_gpickle('./tmp/wtmp'+str(stepCountTmp))
                               self.__lattice=nx.read_gpickle('./tmp/ltmp'+str(stepCountTmp))
                               infile = open('./tmp/WPAtmp'+str(stepCountTmp), 'rb')
                               self.__WalkPosArray[0]=nx.read_gpickle(infile)
                               infile.close()
                               infile = open('./tmp/optmp'+str(stepCountTmp), 'rb')
                               oldPos=nx.read_gpickle(infile)
                               infile.close()
                               infile = open('./tmp/nptmp'+str(stepCountTmp), 'rb')
                               newPos=nx.read_gpickle(infile)
                               infile.close()
                               infile = open('./tmp/nwtmp'+str(stepCountTmp), 'rb')
                               newWalk=nx.read_gpickle(infile)
                               infile.close()
                               infile = open('./tmp/dtotmp'+str(stepCountTmp), 'rb')
                               self.__distancesToOrigin=nx.read_gpickle(infile)
                               infile.close()
                               infile = open('./tmp/satmp'+str(stepCountTmp), 'rb')
                               self.__SA=nx.read_gpickle(infile)
                               infile.close()
                               infile = open('./tmp/saptmp'+str(stepCountTmp), 'rb')
                               self.__SAPath=nx.read_gpickle(infile)
                               infile.close()
                               infile = open('./tmp/ctmp'+str(stepCountTmp), 'rb')
                               self.__Contr=nx.read_gpickle(infile)
                               infile.close()
                               restartCount=+1
                               restartActive=True
                               stepCountMax=stepCountTmp
                   elif stop==True: 
                       stepPossible = True
                       restartActive=True
                             
               if len(self.__distancesToOrigin) < numberOfSteps+1: stepPossible = True
               else: stepPossible = False
               
               if restartActive == False: 
                   #extend the current walk stats
                   self.__WalkPosArray[0] = [oldPos, newPos, newWalk]    
                   d = np.linalg.norm(np.array(newPos,dtype=float) - np.array(seedPos,dtype=float), self._norm)
                   d = d * d
                   self.__distancesToOrigin.append(d)
                   self.__SA.append(self.__saCount)
                   self.__SAPath.append(self.__saPathCount)
                   self.__Contr.append(self.__ContrCount)
                   self.__WalkPosArray[0] = [oldPos, newPos, newWalk]
                       
               # store overall statistics
               self.__allDistancesToOrigin[n] = self.__distancesToOrigin
               self.__allValidNeighbors[n] = self.__validNeighbors
               self.__allContractNeighbors[n] = self.__contractNeighbors
               self.__allSA[n] = self.__SA
               self.__allSAPaths[n] = self.__SAPath
               self.__allContr[n] = self.__Contr

            # extend stats arrays
            self.__allWalksSA.append(self.__allSA)
            self.__allWalksSAPaths.append(self.__allSAPaths)
            self.__allWalksContr.append(self.__allContr)
            self.__allWalksDist.append(self.__allDistancesToOrigin)
            
    

