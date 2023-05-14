####### FindTopo.py - Started 26/11/2020
####### A Python script for finding defects from simulation input
####### Based off the FindDefects.py script to compute topological charge fields instead

from locale import locale_encoding_alias
from matplotlib.pyplot import sca
import numpy as np
import sys
import math

from numpy.compat import py3k
from numpy.lib.scimath import sqrt

### Handle Script Arguments ##############################################################
## Example usage:
# python3 ./FindDefects.py /path/to/frame/input.txt ./defectOut.txt 1 10
print( "Arguments:" )
for arg in sys.argv:
	print( "\t" + arg )
print()
aDATANAME = sys.argv[1] # path to data .txt file
aOUTNAME = sys.argv[2] # path to write dumped defects to
aLAYERHEIGHT = float(sys.argv[3]) # the layer height, acts as a characteristic scale for the script
aGRIDSPACING = int(sys.argv[4]) # how many "slices" to make in the x and y direction for the grid. should be comparible length scale to RodShapedBacterium length

##########################################################################################

### 'Tuning' Constants ###################################################################
EPSILON_LAYER = 0.5 # tolerance (as a proportion of layer height) that you are allowed to be out from a layer by
N_LOCALPARTS = 4 # number of local particles we must have at minimum to consider computing local Q
#TODO: tune these!


##########################################################################################

### Helper Classes #######################################################################

#class for a particle (to make memory management easier)
class Particle:
    def __init__(self, id, length, diam, pX, pY, pZ, oX, oY, oZ):
        # normal params
        self.ID = id
        self.Length = length
        self.Diameter = diam
        self.posX = pX
        self.posY = pY
        self.posZ = pZ
        self.oriX = oX # ori is the unit director in 3D
        self.oriY = oY
        self.oriZ = oZ

    def getArea(self): #compute the 2D area taken by this particle
        return ((self.Length - self.Diameter) * self.Diameter + math.pi * (self.Diameter/2)**2)

    def getVolume(self): #compute the 3D volume taken by this particle
        return ((math.pi * (self.Diameter/2)**2) * (self.Length - self.Diameter) + 4 * (math.pi (self.Diameter/2)**3) / 3)

    def dir2(self): #returns the normalised 2D director
        xyLength = sqrt(self.oriX**2 + self.oriY**2) #length of the 2D director
        return [self.oriX / xyLength, self.oriY / xyLength]

#class for a "node", nodes being what we search for defects on
class Node:
    def __init__(self, pX, pY, pZ):
        self.posX = pX
        self.posY = pY
        self.posZ = pZ
        self.localParts = [] # a list of particles local (but not necessarily neighbours to!) the node

        #topo related bits
        self.charge = None # topological charge
        self.computedQ = False # flag on whether we've computed Q or not
        self.avDir = [None, None] # averaged director based on local particles
        self.Q = [  [None, None],
                    [None, None]] 

    def assignParticle(self, part: Particle):
        self.localParts.append(part)

    def computeDist2(self, part : Particle): # computes the 2D distance to a particle from this node (ignores z!!)
        dx = part.posX - self.posX
        dy = part.posY - self.posY
        return sqrt(dx*dx + dy*dy)

    def computeDist3(self, part : Particle): # computes the 3D distance to a particle from this node
        dx = part.posX - self.posX
        dy = part.posY - self.posY
        dz = part.posZ - self.posZ
        return sqrt(dx*dx + dy*dy + dz*dz)

    def computeLocalQ(self): #computes the local Q tensor by averaging over local directors
        #return NOW if we don't have enough particles
        if len(self.localParts) < N_LOCALPARTS:
            return

        # compute director by averaging local directors
        xDirs = []
        yDirs = []
        for part in self.localParts:
            partDir = part.dir2()
            xDirs.append(partDir[0])
            yDirs.append(partDir[1])
        self.avDir[0] = sum(xDirs) / len(xDirs)
        self.avDir[1] = sum(yDirs) / len(yDirs)

        self.computeQ()

    def evalTopo(self, grid, x: int, y: int): #checks the topology by using neighbouring nodes
        #TODO: write this
        pass

        ## FLOW
        #
        ## check to make sure none of neighbours are "out of range"
            ## return None if so
        ## check to make sure all neighbours have had Q computed
            ## return None if not
        ## construct list containing an ordered loop of neighbours and call topoChargeLocal

    #computes the local topological charge if given a list of particles that serves as a loop (unordered)
    def topoChargeLocal(self, sortedLoop) -> float:
        ##TODO: this needs rewriting to work on nodes instead of particles while keeping sorted loop

        phi = 0.0 #the running total we're using to compute the topo
        #now go through the loop and compute the smallest topological angle
        for i, part in enumerate(sortedLoop):
            nextPart = None #the next particle along in the loop
            if i+1 == len(sortedLoop): #wrap around list if necessary
                nextPart = sortedLoop[0]
            else:
                nextPart = sortedLoop[i+1]

            # from here replicating the MPCD algorithm
            dir1 = part.dir2()
            dir2 = nextPart.dir2()
            dot = dir1[0]*dir2[0]+dir1[1]*dir2[1]
            det = dir1[0]*dir2[1] - dir1[1]*dir2[0]
            ang = math.atan2(abs(det), dot)

            if ang > 0.5*math.pi: #flip as necessary
                dir2[0] = -dir1[0]
                dir2[1] = -dir2[1]
                dot = dir1[0]*dir2[0]+dir1[1]*dir2[1]
                det = dir1[0]*dir2[1] - dir1[1]*dir2[0]
            sign = 1 #get sign
            if det < 0:
                sign = -1
            
            phi += sign*math.atan2(abs(det), dot) #add the angle to the topo counter
        
        return phi # return the topo counter
    
    def computeQ(self): # computes the local Q tensor based of the average director
        #TODO: FINISH this
        for i in range(2):
            for j in range(2):
                self.Q[i][j] = self.avDir[i]*self.avDir[j]
                if i==j: #substract I term
                    self.q[i][j] -= 0.5

                #TODO: scale everything with order parameter S
        self.computedQ = True

    def topoAngleLocal(self, sortedLoop):
        #TODO: code this
        pass

    def dumpNodeData(self, file):
        file.write(str(self.posX) + "\t" + str(self.posY) + "\t" + str(self.posZ) + "\t" + str(self.charge) + "\n")

##########################################################################################

### Core Function ##########################################################################

def singleFrame(inPath : str, outPath : str, step : int, outType : str):
    #construct infile name
    inFilePath = inPath + str(step).zfill(5) + "." + outType

    #load data file into memory and store as a particle list
    print("Reading file " + inFilePath + " for data.")
    particleList = []
    maxZ = 0 #var we need for later, when doing layers
    try:
        data = np.genfromtxt(inFilePath, skip_header=1) #get data

        #now loop through and populate particle list
        for row in data:
            if row[5] > maxZ:
                maxZ = row[5]
            particleList.append(Particle(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8]))
    except: #throw error and terminate if we can't read for some reason
        print("Failed to read data, terminating.")
        quit()
    totParticleNo = len(particleList)
    #print("Finished loading data. Total particle number: " + str(totParticleNo) + "\n")

    #sort into layers
    layerContainer = [] #the object we'll be using to store our lists of particles
    layerNo = round(maxZ / aLAYERHEIGHT) # number of layers we expect to use
    for i in range(layerNo): # set up each layer's list
        layerContainer.append([])
    #print("Splitting into layers. Expecting " + str(layerNo) + " layers.")

    #loop through particle list and place them into appropriate list based on 
    for part in particleList:
        scaledHeight = part.posZ / aLAYERHEIGHT
        targetLayer = int(scaledHeight) # the integer part only

        if abs(scaledHeight-targetLayer) < EPSILON_LAYER: #if you're within epsilon of the target layer
            layerContainer[targetLayer].append(part) #move particle to the appropriate layer

    # print("Split into the following distribution of layers:")
    # totLayerAssigned = 0 # running total of how many particles ahve been assigned
    # for i in range(len(layerContainer)):
    #     currLayer = len(layerContainer[i])
    #     totLayerAssigned += currLayer #add current layer count to running total
    #     print("\tLayer " + str(i) + " - " + str(currLayer) + " particles")
    # print("\tUnassigned: " + str(totParticleNo - totLayerAssigned))

    #now apply the algorithm for each layer
    for currLayer, layer in enumerate(layerContainer):
        #print("\nComputing defects on layer " + str(currLayer))

        #first, need to find the min and max on x and y coords to construct a bounding box
        maxX, maxY, minX, minY = 0, 0, 0, 0
        for part in layer:
            if part.posX < minX:
                minX = part.posX
            if part.posX > maxX:
                maxX = part.posX
            if part.posY < minY:
                minY = part.posY
            if part.posY > maxY:
                maxY = part.posY
        
        #using this we can construct our grid parameters
        #grid is assumed to start at (minX, minY) and has a new element every dX, dY units
        xGridSize = math.ceil((maxX - minX)/aGRIDSPACING)
        yGridSize = math.ceil((maxY - minY)/aGRIDSPACING)
        grid = [] # init our grid object
        for i in range(xGridSize+1): # aXGRIDSIZE+1 because we want to consider the far edge too
            row = []
            for j in range(yGridSize+1):
                row.append(Node(minX + aGRIDSPACING*i, minY + aGRIDSPACING*j, currLayer * aLAYERHEIGHT))  #construct our grid
            grid.append(row)
        
        #do some binning to assign particles to nearest nodes
        for part in layer:
            # get which cell this particle should be assigned to (closest node)
            cellX = round((part.posX - minX)/aGRIDSPACING) 
            cellY = round((part.posY - minY)/aGRIDSPACING)
            grid[cellX][cellY].assignParticle(part)

        #now compute local Qs
        for y, row in enumerate(grid):
            for x, node in enumerate(row):
                node.computeLocalQ()

        #from local Qs, now compute topo charge field
        topoField = []
        for y, row in enumerate(grid):
            rowTopo = []
            for x, node in enumerate(row):
                rowTopo.append(node.evalTopo(grid, x, y))
            topoField.append(rowTopo)
        
        #handle topo file writing
        outFilePath = outPath + str(step).zfill(5) + "." + outType #construct outfile path
        open(outFilePath, "w").close() #force clear the file
        with open(outFilePath, "w") as file:
            #prepare file with headers
            file.write("posX\tposY\tposZ\tcharge\n")
            #now loop through and dump defects
            for row in grid:
                for node in row:
                    node.dumpNodeData(file)

        #TODO: render and save topo field data for this timestep

singleFrame("./data/vis_biofilm_", "./data/topoout_", 560, "txt")
singleFrame("./data/vis_biofilm_", "./data/topoout_", 580, "txt")
singleFrame("./data/vis_biofilm_", "./data/topoout_", 601, "txt")