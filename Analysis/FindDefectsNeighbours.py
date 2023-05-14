####### FindDefects.py - Started 24/11/2020
####### A Python script for finding defects from simulation input

from locale import locale_encoding_alias
from matplotlib.pyplot import sca
import numpy as np
import sys
import math
import copy

from numpy.compat import py3k
from numpy.lib.scimath import sqrt

### Helper Classes #######################################################################

#class for a particle (to make memory management easier)
# Rory: updated the area and volume calculations - the length is already pole to pole
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
        self.nList = []

    def getArea(self): #compute the 2D area taken by this particle
        return (self.Length * self.Diameter + math.pi * (self.Diameter/2)**2)

    def getVolume(self): #compute the 3D volume taken by this particle
        return ((math.pi * (self.Diameter/2)**2) * self.Length + 4 * (math.pi * (self.Diameter/2)**3) / 3)

    def dir2(self): #returns the normalised 2D director
        xyLength = sqrt(self.oriX**2 + self.oriY**2) #length of the 2D director
        return [self.oriX / xyLength, self.oriY / xyLength]

    def getQ():
        #TODO: code this
        pass

class PointParticle(Particle):
     def __init__(self, id, diam, pX, pY, pZ, oX, oY, oZ):
        super().__init__(id, 0, diam, pX, pY, pZ, oX, oY, oZ)

#class for a defect
class Defect:
    def __init__(self, pX, pY, pZ, charge, ori):
        self.posX = pX
        self.posY = pY
        self.posZ = pZ
        self.charge = charge
        self.orientation = ori

    def getDefectData(self):
        return [self.posX,
                self.posY,
                self.posZ,
                self.charge,
                self.orientation]

    def dump(self, id : int, file):
        file.write(str(id) + "\t" + str(self.posX) + "\t" + str(self.posY) + "\t" + str(self.posZ) + "\t" + str(self.charge) + "\t" + str(self.orientation) + "\n")

#class for a "node", nodes being what we search for defects on
class Node:

    def __init__(self, pX, pY, pZ):
        self.posX = pX
        self.posY = pY
        self.posZ = pZ
        self.nList = [] # neighbours are assumed to be given to the node later
        self.localParts = [] # a list of particles local (but not necessarily neighbours to!) the node
        self.localIDs = {} #dictionary containing local IDs
        self.localDensity = 0 # local density we need to compute

    def computeNeighbours(self, pList):
        for particle in pList:
            dist = self.computeDist2(particle)
            if dist < R_NEIGHBOUR*aLAYERHEIGHT: # if you're in the "neighbourhood" region, then add to nList
                self.nList.append(particle)

    # set the "target" particle (part) as a local particle if it isn't in the list of local particles already
    def setLocalPart(self, part: Particle):
        ID = part.ID
        if ID in self.localIDs:
            return # if the target ID is already in our dictionary then don't add it again

        self.localIDs[ID] =  self.localParts.append(part) #if we're here then the target isn't in the list of local parts, so add it

    # assign the particle, and it's neighbours, as a local particle
    def assignParticle(self, part: Particle):
        self.setLocalPart(part)

        for p in part.nList: #try add all neighbours of this particle to the local list
            self.setLocalPart(p)

    def assignNeighbours(self):
        self.computeNeighbours(self.localParts) # compute neighbours from the local cell list
        self.computeDensity()

    def computeDensity(self):
        runningTotalDens = 0 #running total of the density as we calculate it
        for particle in self.nList: # loop through particles
            dist = self.computeDist2(particle) # compute the distance
            if dist <= R_DENSITY*aLAYERHEIGHT: # if you're in the region we're considering of density, note the area
                runningTotalDens += particle.getArea()

        # set the density as being the proportion of the running total to the circle made by the area we're considering
        self.localDensity = runningTotalDens / (math.pi * (R_DENSITY**2))

    def computeDist2(self, part : Particle): # computes the 2D distance to a particle from this node (ignores z!!)
        dx = part.posX - self.posX
        dy = part.posY - self.posY
        return sqrt(dx*dx + dy*dy)

    def computeDist3(self, part : Particle): # computes the 3D distance to a particle from this node
        dx = part.posX - self.posX
        dy = part.posY - self.posY
        dz = part.posZ - self.posZ
        return sqrt(dx*dx + dy*dy + dz*dz)

    def evalTopo(self): #checks the topology of this node
        #construct lengths to each particle
        dist = []
        for part in self.localParts:
            dist.append(self.computeDist2(part))

        for R in np.arange(R_NEIGHBOUR, aSEARCHSPACING, -aSEARCHSPACING):
            loop = self.constructLoop(R, dist) #construct a loop
            if loop != 0: #if we've constructed a loop with what we think is enough elements
                loop = self.orderLoop(loop) #sort the loop
                charge = self.topoChargeLocal(loop) #now try computing the charge
                if abs(abs(charge)-0.5) < EPSILON_TOPOTHRESHOLD: # if you have \pm 1/2 charge
                    angle = self.topoAngleLocal(loop) #compute angle
                    return Defect(self.posX, self.posY, self.posZ, charge, angle)

        return 0 # if we're here then we clearly don't have a defect, dump 0 to indicate this

    #search for a loop of sufficient size for a given radius R
    def constructLoop(self, R : float, dist):
        potentialLoop = [] # construct our potential loop
        for i, part in enumerate(self.localParts):
            if dist[i] > R - epsLoop and dist[i] < R + epsLoop:
                potentialLoop.append(part)

        #now check if the potential loop is sufficient, return as necessary
        if len(potentialLoop) >= N_LOOPMIN:
            return potentialLoop

        return 0 #if we get here we've found no sufficient loop

    def orderLoop(self, unsortedLoop):
        #need to sort this list into an ordered loop, do so by using the angle to each particle
        angles = [] #compute the angle each particle makes with with the node
        for part in unsortedLoop:
            dx = part.posX - self.posX
            dy = part.posY - self.posY
            angles.append(math.atan2(dy, dx))
        #now sort the loop by angles
        return [x for _, x in sorted(zip(angles, unsortedLoop))]

    #computes the local topological charge if given a list of particles that serves as a loop (unordered)
    def topoChargeLocal(self, sortedLoop) -> float:
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
                dir2[0] = -dir2[0]
                dir2[1] = -dir2[1]
                dot = dir1[0]*dir2[0]+dir1[1]*dir2[1]
                det = dir1[0]*dir2[1] - dir1[1]*dir2[0]
            sign = 1 #get sign
            if det < 0:
                sign = -1

            phi += sign*math.atan2(abs(det), dot) #add the angle to the topo counter

        return 0.5*phi/math.pi # return the topo counter

    def topoAngleLocal(self, sortedLoop):
        #TODO: code this
        pass

##########################################################################################

### Core Function ##########################################################################


# alt. core fn to interact directly with data loaded in visualise_biofilm.py
# this function will assume particleList is already filtered for the relevent layer
def singleFrameFromCellData(cellList, gridSpacing : float, layerZ : float):

    particleList=[]
    for cell in cellList:
        artificialPoints, diam = genAdditionalFromCell(cell)
        for artificialPoint in artificialPoints:
            # make cell into a uniformly populated point cloud
            # each point inside the cell is given the same orientaion as the original
            particleList.append(
                PointParticle(cell.cell_id,diam,*artificialPoint,*cell.ori)
                )

    #first, need to find the min and max on x and y coords to construct a bounding box
    maxX, maxY, minX, minY = 0, 0, 0, 0
    for part in particleList:
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
    xGridSize = math.ceil((maxX - minX)/gridSpacing)
    yGridSize = math.ceil((maxY - minY)/gridSpacing)
    grid = [] # init our grid object
    for i in range(xGridSize+1): # aXGRIDSIZE+1 because we want to consider the far edge too
        row = []
        for j in range(yGridSize+1):
            row.append(Node(minX + gridSpacing*i, minY + gridSpacing*j,layerZ))  #construct our grid
        grid.append(row)

    #do some binning to assign particles to nearest nodes
    for part in particleList:
        # get which cell this particle should be assigned to (closest node)
        cellX = round((part.posX - minX)/gridSpacing)
        cellY = round((part.posY - minY)/gridSpacing)
        grid[cellX][cellY].assignParticle(part)

    #now perform the actual computation
    for y, row in enumerate(grid):
        for x, node in enumerate(row):
            #first compute neighbours from the particles binned in the last step
            node.assignNeighbours()

            #only consider computing defects if we are abowve a critical density in this node
            if (node.localDensity > DELTA_DENSITY):
                currDefect = node.evalTopo()
                if currDefect != 0: # if we didn't get the value 0 back, then we have a defect!
                    defectContainer.append(currDefect)

    # Can we evaluate topological charge here too?
    # #now compute local Qs
    # for y, row in enumerate(grid):
    # 	for x, node in enumerate(row):
    # 		node.computeLocalQ()
    #
    # #from local Qs, now compute topo charge field
    # topoField = []
    # for y, row in enumerate(grid):
    # 	rowTopo = []
    # 	for x, node in enumerate(row):
    # 		rowTopo.append(node.evalTopo(grid, x, y))
    # 	topoField.append(rowTopo)
    # return np.asarray([[ node.getChargeData() for node in row ] for row in grid])
    return np.asarray([defect.getDefectData() for defect in defectContainer])

def singleFrame(inPath : str, neighbourInPath : str, outPath : str, step : int, outType : str):
    #construct infile name
    inFilePath = inPath + str(step).zfill(5) + "." + outType
    neighbourInFilePath = neighbourInPath + str(step).zfill(5) + "." + outType

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
    layerNo = int(round(maxZ / aLAYERHEIGHT)) # number of layers we expect to use
    for i in range(layerNo): # set up each layer's list
        layerContainer.append([])
    #print("Splitting into layers. Expecting " + str(layerNo) + " layers.")

    #load in neighbour data and link to each corresponding particle
    print("Reading file " + neighbourInFilePath + " for neighbour data.")
    try:
        #need to manually construct list because of varying column lengths, thanks numpy
        data = [] #empty data structure
        f = open(neighbourInFilePath, 'r')
        while f:
            line = f.readline().rstrip("\t").rstrip("\n")
            if(not line):
                break
            row = line.split("\t")
            data.append(row)
        f.close()

        #construct hash list
        lst_to_dist = lambda lst: {
            int(lst[ii].ID): lst[ii] for ii in range(0,len(lst))
        }
        cell_hash = copy.copy(lst_to_dist(particleList))

        print("Neighbour ID dictionary constructed, assigning global neighbours")
        #now, go through each row in the data and append neighbours
        for i, row in enumerate(data):
            for j,_ in enumerate(row):
                cell_hash[int(row[0])].nList.append(cell_hash[int(row[j])])
    except: #throw error and terminate if we can't read for some reason
        print("Failed to read neighbour data, terminating.")
        quit()
    print("Finished reading neighbour data")

    #loop through particle list and place them into appropriate list based on
    for part in particleList:
        scaledHeight = part.posZ / aLAYERHEIGHT
        targetLayer = int(scaledHeight) # the integer part only

        if abs(scaledHeight-targetLayer) < EPSILON_LAYER: #if you're within epsilon of the target layer
            try:
                layerContainer[targetLayer].append(part) #move particle to the appropriate layer
            except:
                pass

    #now want to construct neighbours WITHIN EACH LAYER ONLY - If a particle isn't in a the same layer, despite being neighbours, then don't add them to the same nList
    #first, load the neighbour data

    # print("Split into the following distribution of layers:")
    # totLayerAssigned = 0 # running total of how many particles ahve been assigned
    # for i in range(len(layerContainer)):
    #     currLayer = len(layerContainer[i])
    #     totLayerAssigned += currLayer #add current layer count to running total
    #     print("\tLayer " + str(i) + " - " + str(currLayer) + " particles")
    # print("\tUnassigned: " + str(totParticleNo - totLayerAssigned))


    #now apply the algorithm for each layer
    defectContainer = [] #init the container used for defects
    for layer in layerContainer:
        currLayer = layerContainer.index(layer)
        print("\nComputing defects on layer " + str(currLayer))

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

        #remove particles not in this layer from the appropriate neighbour list
        lst_to_dist = lambda lst: {
            int(lst[ii].ID): lst[ii] for ii in range(0,len(lst))
        }
        cell_hash = copy.copy(lst_to_dist(layer)) #construct a hash table for this layer
        for part in layer:
            for neighbour in part.nList: # for each neighbour...
                if not (neighbour.ID in cell_hash): # if it's NOT in the hash list (ie not in layer)
                    part.nList.remove(neighbour) #remove it

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
            cellX = int(round((part.posX - minX)/aGRIDSPACING))
            cellY = int(round((part.posY - minY)/aGRIDSPACING))
            grid[cellX][cellY].assignParticle(part)

        #now perform the actual computation
        for y, row in enumerate(grid):
            for x, node in enumerate(row):
                #first compute neighbours from the particles binned in the last step
                #FIXME: perform the calculation necessary to get around the above fixme
                node.assignNeighbours()

                #only consider computing defects if we are abowve a critical density in this node
                if (node.localDensity > DELTA_DENSITY):
                    currDefect = node.evalTopo()
                    if currDefect != 0: # if we didn't get the value 0 back, then we have a defect!
                        defectContainer.append(currDefect)

    #handle file writing
    outFilePath = outPath + str(step).zfill(5) + "." + outType #construct outfile path
    open(outFilePath, "w").close() #force clear the file
    with open(outFilePath, "w") as file:
        #prepare file with headers
        file.write("id\tposX\tposY\tposZ\tcharge\torientation\n")
        #now loop through and dump defects
        if len(defectContainer) > 0:
            for i, defect in enumerate(defectContainer):
                defect.dump(i, file)
            print("\nDumped defects to " + outFilePath + "\n")
        else:
            print("\nNo defects found for step " + str(step) + "!")

#### Program flow notes ####
#
## Read in data
## Sort data into layers, for each layer...
    ## find the min and max on x and y coords, to serve as the bounding box
    ## construct a grid based on the above
    ## compute neighbour lists on nodes within layer
        ## use this to compute a local density value (just sum width * lengths of neighbours)
    ## for each particle, compute order param and other LC stuff
    ## for each node, if density is sufficient....
        ## compute "loops" using neighbour lits (either try a recursive node based search or do the dumb epsilon method)
            ## go around each "loop" to compute the topo
                ## increase loop width if no topo found up to a certain limit (likely a function of density and average length of bacteria in this layer)
                ## if topo found then compute topo orientation
## dump topo + orientation list


if __name__ == "__main__":
    ### Handle Script Arguments ##############################################################
    ## Example usage:
    # python3 ./FindDefects.py /path/to/frame/input.txt /path/to/neighbout/input.txt ./defectOut.txt 1 10 2
    print( "Arguments:" )
    for arg in sys.argv:
        print( "\t" + arg )
    print()
    aDATANAME = sys.argv[1] # path to data .txt file
    aNEIGHBOURDATANAME = sys.argv[2] # path to neighbour data .txt file
    aOUTNAME = sys.argv[3] # path to write dumped defects to
    aLAYERHEIGHT = float(sys.argv[4]) # the layer height, acts as a characteristic scale for the script
    aGRIDSPACING = int(sys.argv[5]) # how many "slices" to make in the x and y direction for the grid. should be comparible length scale to RodShapedBacterium length
    aSEARCHSPACING = float (sys.argv[6]) # how much to vary spacing when searching for defects

    ##########################################################################################

    ### 'Tuning' Constants ###################################################################
    EPSILON_LAYER = 0.1 # tolerance (as a proportion of layer height) that you are allowed to be out from a layer by
    EPSILON_LOOP = 0.05 # tolerance (as a proportion of grid spacing) that you are allowed to be out from a loop by
    EPSILON_TOPOTHRESHOLD = 0.00001 # tolerance of how far you can be out from \pm 1/2 to be considered having half topological charge
    R_NEIGHBOUR = 20 # how many units of aLAYERHEIGHT to consider out to when computing neighbour lists.
        # This will also constrain the size of the toposearch!!!
    R_DENSITY = 20 # how many units of aLAYERHEIGHT to consider when computing local densities. Should be less than R_NEIGHBOUR
    DELTA_DENSITY = 0.1 # only try computing a defect if node density is higher than this
    #TODO: tune these!
    N_LOOPMIN = 25 # recquire a minimum of this amount of particles in a loop to consider it a valid loop

    if(R_DENSITY > R_NEIGHBOUR):
        print("Constant warning! \nR_DENSITY is constrainted by R_NEIGHBOUR. Adjust hardcoded 'tuning' constants!")


    epsLoop = EPSILON_LOOP*aGRIDSPACING

    # Populate Node relevant globals
    Node.aDATANAME = aDATANAME
    Node.aOUTNAME = aOUTNAME
    Node.aLAYERHEIGHT = aLAYERHEIGHT
    Node.aGRIDSPACING = aGRIDSPACING
    Node.aSEARCHSPACING = aSEARCHSPACING
    Node.EPSILON_LAYER = EPSILON_LAYER
    Node.EPSILON_LOOP = EPSILON_LOOP
    Node.EPSILON_TOPOTHRESHOLD = EPSILON_TOPOTHRESHOLD
    Node.R_NEIGHBOUR = R_NEIGHBOUR
    Node.R_DENSITY = R_DENSITY
    Node.DELTA_DENSITY = DELTA_DENSITY
    Node.N_LOOPMIN = N_LOOPMIN
    epsLoop = EPSILON_LOOP*aGRIDSPACING
    Node.epsLoop = epsLoop

    ## Debugging lines #######################################################################
    singleFrame("./data/vis_biofilm_", "./data/neighbours_", "./data/defectout_", 400, "txt")
    #singleFrame("./data/vis_biofilm_", "./data/neighbours_", "./data/defectout_", 419, "txt")
    #singleFrame("./data/vis_biofilm_", "./data/neighbours_", "./data/defectout_", 420, "txt")
