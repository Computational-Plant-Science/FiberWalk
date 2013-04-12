'''
-------------------------------------------------------------------------------------------
FIBER WALK DEMO
This is the main.py file, which demonstrates the use of the FiberWalk simulation class.
It contains basic examples to run the simulation, set the simulation parameters 
and access the data to calculate statistics and a basic plotting example

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

import FiberWalk as fw
import os
import pickle
import networkx as nx
import pylab as p
import numpy as np
from scipy import polyfit, polyval, sqrt

# simple helper function to generate names to save
def saveName():
    fname = '_Contraction' + str(numberOfContractions) + '_Steps' + str(numberOfSteps)
    fpath = './FiberWalks/' + str(dimension) + 'D/'
    print 'Checking directory structure: ' + fpath
    try:
        os.stat(fpath)
    except:
        os.makedirs(fpath)
    return fname, fpath



# Log Log Plot of the distances to origin
def makeAvgDistanceFigure(fname, fpath, allDists, numberOfSteps):
    maxLength = 0
    maxIdx = 0
    minL=numberOfSteps-1
    for i in range(len(allDists)):
        if maxLength < len(allDists[i]):
            maxLength = len(allDists[i])
            maxIdx = i
    
    avgDist = np.zeros(minL, dtype=float)
    numberOfWalksArray = np.zeros(minL, dtype=float)
    
    
    for i in range(len(allDists)):
        if (len(allDists[i]) - 1) >= minL:
            for j in range(minL):
                avgDist[j] = avgDist[j] + allDists[i][j]
                numberOfWalksArray[j] += 1
    
    avgDist = np.array(avgDist) / np.array(numberOfWalksArray)
    allY = []
    for i in range(len(avgDist)):
        allY.append(i)
    logDists = np.log(avgDist[1:])
    logY = np.log(allY[1:])
    (ar, br) = polyfit(logY[30:minL], logDists[30:minL], 1)
    xr = polyval([ar, br], logY[30:minL])
    xr = np.exp(xr)
    #compute the mean square error
    err = sqrt(sum((xr - logDists[30:minL]) ** 2) / len(logDists[30:minL]))
    print('Linear regression using polyfit')
    print('Regression: a=%.2f b=%.2f' % (ar, br))
    
    p.figure(1)
    p.plot(allY[1:minL], avgDist[1:minL], c='b', marker='.', markersize=3)
     
    p.plot(allY[31:minL], xr, c='r',linewidth=2)
    p.ylim([1, 10000])
    p.xlim([1, 10000])
    p.grid(True, which='both')
    p.xscale('log')
    p.yscale('log')
    p.xlabel('time (steps along edges)')
    p.ylabel('MSD from the seed point')
    p.title('Mean Square Distance per step')
    p.text(60, 20, 'y=' + str(np.around(ar, 2)) + 'x+' + str(np.around(br, 2)))
    p.text(50, 30, '#walks: ' + str(numberOfWalksArray[0]))
    print 'save plot to '+fpath
    p.savefig(fpath + fname + '_avg')
    p.close()
    
    

  
#Fiber Walk parameters    
dimension = 2 # choose dimension
numberOfSteps = 100 # choose the length of the walk
numberOfObjects = 5 # choose number of walks to be generated
numberOfContractions = 1 # choose number of contractions per step


print 'Initializing Fiber Walk'
FW = fw.FiberWalk(dimension, numberOfObjects, numberOfContractions)
print 'The Fiber walks ...'
FW.FiberWalk(numberOfSteps)
fname, fpath = saveName()

#getting data
walk=FW.get_walk()
lat=FW.get_lattice()
dists=FW.get_all_distances_to_origin()

# saving data as pickles
nx.write_gpickle(walk, './walk'+fname)
nx.write_gpickle(lat, './lattice'+fname)
outfile = open('./dists'+fname, 'wb')
pickle.dump(dists,outfile)
outfile.close()

# loading data
infile = open('./dists'+fname, 'rb')
distsLoad=nx.read_gpickle(infile)
infile.close()

#show data
print 'Example of showing loaded data'
print distsLoad
#make a plot
makeAvgDistanceFigure(fname, fpath, distsLoad, numberOfSteps)
print 'The Fiber walked ;-)'