"""
    Find defects from simulation input

    Utilities include both finding defects given an input file or dynamically
    given a list of appropriately filtered particles.

    This script is intended for use with any kind particle
    TODO: move the general particle definition and use this to inherit from in
    RodShapedBacterium and ChainingRodShapedBacterium scripts

    Program flow notes
    Read in data
    Sort data into layers, for each layer...
        find the min and max on x and y coords, to serve as the bounding box
        construct a grid based on the above
        compute neighbour lists on nodes within layer
            use this to compute a local density value (just sum width * lengths of neighbours)
        for each particle, compute order param and other LC stuff
        for each node, if density is sufficient....
            compute "loops" using neighbour lits (either try a recursive node based search or do the dumb epsilon method)
                go around each "loop" to compute the topo
                    increase loop width if no topo found up to a certain limit (likely a function of density and average length of bacteria in this layer)
                    if topo found then compute topo orientation
    dump topo + orientation list

    Example usage:
       python3 ./FindDefects.py /path/to/frame/input.txt ./defectOut.txt 1 10 0.5

    Author: Tim Kozhukhov
        FindDefects.py - Started 24/11/2020
    Edited by: Rory Claydon
"""

from locale import locale_encoding_alias
from matplotlib.pyplot import sca
from numpy.compat import py3k
from numpy.lib.scimath import sqrt

# Standard libraries
import sys
import math
import numpy as np
import numpy.linalg
import pandas as pd
from typing import TypeVar, Type
import re
import numpy.random
import matplotlib.pyplot as plt

# User defined
from RodShapedBacteria import RodShapedBacterium
from ChainingRodShapedBacteria import ChainingRodShapedBacterium
from AG43RodShapedBacteria import AG43RodShapedBacterium

"""
    Utility functions
"""

# Move this somewhere more sensible later
def intialiseElementsFromData(data,class_dict):
    """
        Populate a list of elements from a pandas data frame
        Parameters:
            data: pandas dataframe
                all required data to initialise members of type Class
            class_dict: dict
                These define how objects will appear on the axes. The list should
                be a dictionary of cell_type keys with value the class to use to
                represent this.
        Returns:
            List of members
    """
    cells = []
    species_load_dict={}
    for index in data.index:
        try:
            class_name=data.loc[index,:].cell_type
            Class=class_dict[class_name]
            try:
                cells.append( Class(**data.iloc[index,1:].to_dict()) )
            except Exception as e:
                # print(e)
                try:
                    param_dict=data.iloc[index,1:7].to_dict()
                    param_dict.pop('length', None)
                    cells.append( Class(**param_dict) )
                except:
                    raise
        except Exception as e:
            # print(e)
            # print(f"Failed loading class {class_name}, defaulting to base class RodShapedBacterium")
            cells.append( RodShapedBacterium(**data.iloc[index,1:10].to_dict()) )

        # Keep track of how many of each type are present in the colony
        if class_name not in species_load_dict:
            species_load_dict[class_name]=1
        else:
            species_load_dict[class_name]+=1

    try:
        print("Created chaining hash list")
        RodShapedBacterium.makeHashList(cells)
        # AG43RodShapedBacterium.createSpringHashList(cells)
    except Exception as e:
        print(e)
        quit()
    return cells, species_load_dict


def cellPlaneIntersect(cell,plane_loc,eps,select):
    """
        Determine if any part of the main axis of the cell intersects with the
        selected plane

        Parameters:
            cell: Cell Class instance
                cell to check for instersection
            plane_loc: float
                the position of the plane in the selected axis
            eps: float
                tolerance within which cells will still be included
            select: int
                the selected axis which is projected onto. Alternatively the
                normal of the plane. x,y,z for 0,1,2 resp.

        Returns:
            interaction: bool
                true if cell intersects plane within tolerance
    """

    # Take length from tip to tip
    full_length = 1 + cell.length

    # Plane normal is in the direction of projection
    # using in case of generalistion later
    plane_normal = np.zeros(3)
    plane_normal[select] = 1

    # find the length of the cell projected along the normal to the plane
    cell_proj_len = full_length*cell.ori.dot(plane_normal)

    # Check if any part intersects with plane and if so return true
    mu_upper = ( plane_loc+eps-cell.rcm.dot(plane_normal) ) / cell_proj_len
    mu_lower = ( plane_loc-eps-cell.rcm.dot(plane_normal) ) / cell_proj_len

    if abs(mu_upper)<=0.5 or abs(mu_lower)<=0.5:
        # Limits intersect cell axis
        return True
    elif mu_upper*mu_lower<0:
        # cell is between the two limits
        return True
    else:
        # No direct intersection or bounding of the cell
        if (abs(cell.rcm[select] - plane_loc) < eps):
            center_top    = cell.rcm + 0.5*full_length*cell.ori
            center_bottom = cell.rcm - 0.5*full_length*cell.ori
            print("instersection lag failed")
            print(f"plane at {plane_loc} cell rcm proj {cell.rcm[select]} {center_bottom[select]}")
            print(f"plane lims {plane_loc-eps} {plane_loc+eps}")
            print(f"{( plane_loc-eps-center_bottom.dot(plane_normal) ) / cell_proj_len}")
            quit()
        return False

"""
    User defined types
"""

Particle = TypeVar("Particle",RodShapedBacterium,ChainingRodShapedBacterium)

"""
    Utility Class Definitions
"""

# class RodParticle:
#     """
#         Representation of spherocylinder particle
#     """
#     def __init__(self, id, length, diam, pX, pY, pZ, oX, oY, oZ):
#         self.id = id            # Particle unbique identifier
#         self.length = length    # pole to pole length of the particle
#         self.diameter = diam    # rod diameter
#         self.posX = pX          # pX - position of the com
#         self.posY = pY
#         self.posZ = pZ
#         self.oriX = oX # ori is the unit director in 3D
#         self.oriY = oY
#         self.oriZ = oZ
#
#         self.rcm = np.array([posX,posY,posZ])
#         self.ori = np.array([oriX,oriY,oriZ])
#
#     def getArea(self):
#         """
#             Compute the area of the projection of this particle onto the specified
#             plane
#         """
#         in_plane_length = self.ori[self.tangent_axes]*self.length
#         in_plane_length = np.linalg.norm(in_plane_length)
#         radius = 0.5*self.diameter
#         return in_plane_length*self.diameter + np.pi*radius**2
#
#     def getVolume(self):
#         """
#             Compute the 3D volume of this particle
#         """
#         radius = 0.5*self.diameter
#         return np.pi*self.length*radius**2 + (4*np.pi*radius**3)/3
#
#     def getInPlaneDirector(self):
#         """
#             Find the normalised in plane director
#         """
#         director2D = self.ori[self.tangent_axes]
#         return getUnitVec(director2D)
#
#     def getQ():
#         #TODO: code this
#         pass
#
# class PointParticle(Particle):
#      def __init__(self, id, diam, pX, pY, pZ, oX, oY, oZ):
#         super().__init__(id, 0, diam, pX, pY, pZ, oX, oY, oZ)

class Defect:
    """
        Representation of defects
    """
    def __init__(self, pX, pY, pZ, charge, ori):
        self.posX = pX
        self.posY = pY
        self.posZ = pZ
        self.charge = charge
        self.orientation = ori # should this be a vector?

        self.pos = np.array([pX,pY,pZ])

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}".format(*self.pos,self.charge,self.orientation)

    def getDefectData(self):
        return [*self.pos,
                self.charge,
                self.orientation]

    def dump(self, id : int, file):
        file.write(f"{id}\t{self.__str__()}\n")
        # file.write(str(id)
        #             + "\t"
        #             + str(self.posX)
        #             + "\t"
        #             + str(self.posY)
        #             + "\t"
        #             + str(self.posZ)
        #             + "\t"
        #             + str(self.charge)
        #             + "\t"
        #             + str(self.orientation)
        #             + "\n")


class Node:
    """
        Grid nodes for reference points. Defects are are searched for using each
        Node's neighbouring particles
    """
    def __init__(self, pX, pY, pZ):
        self.posX = pX  # Node location
        self.posY = pY
        self.posZ = pZ
        self.pos = np.array([pX,pY,pZ])
        self.nList = [] # neighbours are assumed to be given to the node later
        self.local_parts = [] # list of particles local (but not necessarily neighbours to!) the node
        self.localDensity = None # local density we need to compute

    @classmethod
    def setHyperParams(cls,
                       # dataname,
                       # outname,
                       layerheight=1,
                       gridspacing=10,
                       searchspacing=0.5,
                       epsilon_layer=0.3,
                       epsilon_loop=0.1,
                       epsilon_topothreshold=0.03,
                       r_neighbour=7,
                       r_density=5,
                       delta_density=0.25,
                       n_loopmin=5):
        """
            Set Node input/output options and defect computation hyperparameters

            Parameters:
            # dataname: str
            #     input data file name
            # outname: str
            #     output data file name
            layerheight: float
                the layer height (acts as a characteristic scale for the script)
            gridspacing: float
                how many "slices" to make in the x and y direction for the grid.
                Should be comparible length scale to RodShapedBacterium length
            searchspacing: float
                how much to vary spacing when searching for defects (min loop radius)
            epsilon_layer: float
                tolerance (as a proportion of layer height) that you are allowed
                to be out from a layer by
            epsilon_loop: float
                tolerance (as a proportion of grid spacing) that you are allowed
                to be out from a loop by
            epsilon_topothreshold: float
                tolerance of how far you can be out from \pm 1/2 to be considered
                having half topological charge
            r_neighbour: float
                how many units of layerheight to consider out to when computing
                neighbour lists. (max radius)
            r_density: float
                how many units of alayerheight to consider when computing local
                densities. Should be less than r_neighbour
            delta_density: float
                 only try computing a defect if node density is higher than this
            n_loopmin: int
                require a minimum of this amount of particles in a loop to
                consider it a valid loop
        """

        # #----------------------- Set'I/O' Constants ----------------------------
        # cls.dataname = dataname
        # cls.outname = outname

        #----------------------- Set length scales -----------------------------
        cls.layerheight = layerheight
        cls.gridspacing = gridspacing
        cls.searchspacing = searchspacing

        #---------------------- Set'Tuning' Constants --------------------------
        cls.epsilon_layer = epsilon_layer
        cls.epsilon_loop = epsilon_loop
        cls.eps_loop = epsilon_loop*gridspacing
        cls.epsilon_topothreshold = epsilon_topothreshold
        cls.r_neighbour = r_neighbour
        cls.r_density = r_density #
        cls.delta_density = delta_density
        cls.n_loopmin = n_loopmin

        # if(r_neighbour*sqrt(2) > gridspacing):
        #     print(
        #     """Constant warning!
        #        r_neighbour is constrainted by agridspacing.
        #        Adjust hardcoded 'tuning' constants or increase grid spacing!
        #        Exiting...
        #        """)
        #     quit()
        if(r_density > r_neighbour):
            print(
            """Constant warning!
               r_density is constrainted by r_neighbour.
               Adjust hardcoded 'tuning' constants!
               Exiting...
               """)
            quit()

    def computeNeighbours(self, part_list):
        """
            Compute the list of neighbours to this Node.

            Parameters:
                part_list: list of particle instances

            Effect:
                updates attribute nList with the this Node's neighbour particles
        """
        #FIXME: Currently this class is being designed for grid based geometries
        # because of how slow the below function is.
        # Would be nice to speed it up and allow arbitrary Node geometry
        for particle in part_list:
            dist = self.computeDist2D(particle)
            if dist < Node.r_neighbour*Node.layerheight:
                # if you're in the "neighbourhood" region, then add to nList
                self.nList.append(particle)

    def assignParticle(self, part: Particle):
        self.local_parts.append(part)

    def assignNeighbours(self):
        # compute neighbours from the local cell list
        self.computeNeighbours(self.local_parts)
        self.computeDensity()

    def computeDensity(self):
        """
            Find the 2D projected area
        """
        runningTotalDens = 0
        for particle in self.nList: # loop through particles
            dist = self.computeDist2D(particle) # compute the distance
            # if you're in the region we're considering of density, note the area
            if dist <= Node.r_neighbour*Node.layerheight:
                runningTotalDens += particle.getArea()

        # set the density as being the proportion of the running total to the
        # circle made by the area we're considering
        self.localDensity = runningTotalDens / (np.pi*Node.r_density**2)

    def computeDist2D(self, part : Particle):
        """
            update to projected version later
        """
        # dx = part.posX - self.posX
        # dy = part.posY - self.posY
        # return sqrt(dx*dx + dy*dy)
        tangent_axes=[0,1]
        return np.linalg.norm(part.rcm[tangent_axes] - self.pos[tangent_axes])

    def computeDist3D(self, part : Particle):
        # dx = part.posX - self.posX
        # dy = part.posY - self.posY
        # dz = part.posZ - self.posZ
        # return sqrt(dx*dx + dy*dy + dz*dz)
        return np.linalg.norm(part.rcm - self.pos)

    def evalTopo(self):
        """
            checks the topology of this node
        """
        #construct lengths to each particle
        dist = []
        for part in self.local_parts:
            dist.append(self.computeDist2D(part))

        for R in np.arange(Node.searchspacing,Node.r_neighbour,Node.searchspacing):
            loop = self.constructLoop(R, dist) # construct a loop
            if loop != None:

                # if we've constructed a loop with what we think is enough elements
                loop = self.orderLoop(loop) # sort the loop

                # now try computing the charge
                charge = self.topoChargeLocal(loop)

                # if you have \pm 1/2 charge
                if abs(abs(charge)-0.5) < Node.epsilon_topothreshold:
                    angle = self.topoAngleLocal(loop) #compute angle
                    return Defect(*self.pos,charge,angle)

        # No defects found
        return None

    def constructLoop(self, R : float, dist):
        """
            Search for a loop of sufficient size for a given radius R
        """

        potentialLoop = [] # construct our potential loop
        for i, part in enumerate(self.local_parts):
            if dist[i] > R - Node.eps_loop and dist[i] < R + Node.eps_loop:
                potentialLoop.append(part)

        #now check if the potential loop is sufficient, return as necessary
        if len(potentialLoop) >= Node.n_loopmin:
            return potentialLoop
        else:
            return None

    def orderLoop(self, unsorted_loop):
        """
            Sort particles in the contour using the angles their rcm makes with
            the Node i.e. polar coordinate (r,theta) with r=|cell.rcm-Node.pos|
            and theta is arctan2(r_1/r_0) with r_i=(cell.rcm-Node.pos)[i]
        """
        # angles = [] # compute the angle each particle makes with with the node
        # for part in unsorted_loop:
        #     dx = part.posX - self.posX
        #     dy = part.posY - self.posY
        #     angles.append(math.atan2(dy, dx))
        get_r = lambda part: part.rcm - self.pos
        sorted_loop = sorted(unsorted_loop,
                             key=lambda cell: math.atan2(*get_r(cell)[[1,0]])
                             )
        # #now sort the loop by angles
        # return [x for _, x in sorted(zip(angles, unsorted_loop))]
        return sorted_loop

    def topoChargeLocal(self, sorted_loop) -> float:
        """
            computes the local topological charge if given a list of particles
            that serves as a loop (unordered)
        """
        # the running total we're using to compute the topo
        phi = 0.0
        for i, part in enumerate(sorted_loop):
            # compute the smallest topological angle

            # wrap around
            nextPart = sorted_loop[(i+1)%len(sorted_loop)]

            # Find the angle as a function of arc length around the contour
            dir1 = part.getInPlaneDirector()
            dir2 = nextPart.getInPlaneDirector()

            # cos(theta)
            d1_dot_d2 = np.dot(dir1,dir2)

            # get the normal to dir2
            perp_to_d2 = np.array([dir2[1],-dir2[0]])

            # sin(theta)
            det = np.dot(dir1,perp_to_d2)
            ang = math.atan2(abs(det), d1_dot_d2)

            if ang > 0.5*math.pi: #flip as necessary
                # dir2[0] = -dir2[0]
                # dir2[1] = -dir2[1]
                dir2 = -dir2
                # dot = dir1[0]*dir2[0]+dir1[1]*dir2[1]
                # det = dir1[0]*dir2[1]-dir1[1]*dir2[0]

                det = -det
                d1_dot_d2 = -d1_dot_d2

            # add the angle to the topo counter
            phi += np.sign(det)*math.atan2(abs(det),d1_dot_d2)

        return 0.5*phi/math.pi # return the topo counter

    def topoAngleLocal(self, sorted_loop):
        #TODO: code this
        pass


"""
    Core functions for defect calculation
"""
# alt. core fn to interact directly with data loaded in visualise_biofilm.py
# this function will assume particle_list is already filtered for the relevent layer
def singleFrameFromCellData(particle_list,plane_loc,
                            hyperparameter_dict=None):

    # particle_list=[]
    # for cell in cellList:
    #     artificialPoints, diam = genAdditionalFromCell(cell)
    #     for artificialPoint in artificialPoints:
    #         # make cell into a uniformly populated point cloud
    #         # each point inside the cell is given the same orientaion as the original
    #         particle_list.append(
    #             PointParticle(cell.cell_id,diam,*artificialPoint,*cell.ori)
    #             )

    defect_container=[]

    if hyperparameter_dict is not None:
        # User would like to change the hyperparameters
        Node.setHyperParams(**hyperparameter_dict)

    #first, need to find the min and max on x and y coords to construct a bounding box
    maxX, maxY, minX, minY = 0, 0, 0, 0
    for part in particle_list:
        if part.rcm[0] < minX:
            minX = part.rcm[0]
        if part.rcm[0] > maxX:
            maxX = part.rcm[0]
        if part.rcm[1] < minY:
            minY = part.rcm[1]
        if part.rcm[1] > maxY:
            maxY = part.rcm[1]

    #using this we can construct our grid parameters
    #grid is assumed to start at (minX, minY) and has a new element every dX, dY units
    x_gridsize = math.ceil((maxX - minX)/Node.gridspacing)
    y_gridsize = math.ceil((maxY - minY)/Node.gridspacing)
    grid = [] # init our grid object
    for i in range(x_gridsize+1):
        # x_gridsize+1 because we want to consider the far edge too
        row = []
        for j in range(y_gridsize+1):
            # construct our grid - the z height of the grid is always the
            # center of the layer we find the defects for
            row.append(
                Node(minX + Node.gridspacing*i,
                     minY + Node.gridspacing*j,
                     plane_loc)
                     )
        grid.append(row)

    # Bin particles to nearest nodes
    for part in particle_list:
        # get which cell this particle should be assigned to (closest node)
        cellX = int(np.around((part.rcm[0] - minX)/Node.gridspacing))
        cellY = int(np.around((part.rcm[1] - minY)/Node.gridspacing))
        grid[cellX][cellY].assignParticle(part)

    # Just to check nonzero contributions
    # for y, row in enumerate(grid):
    #     for x, node in enumerate(row):
    #         if (len(grid[y][x].local_parts)>0):
    #             print(f"grid[{y}][{x}] has {len(grid[y][x].local_parts)} particles")
        ##FIXME: THERE IS A FAILURE IN LOGIC HERE!!!!!!!!!
            # This means that there will be, at most, a dX*dY square being considered for each node
            # This means that the "scale" of a*GRIDSIZE constrains the neighbour lists!!!
            # Need to actually consider the nearest R_NEIGHBOUR + 1 node's local_parts, not just the individual one

    # now perform the actual computation
    for y, row in enumerate(grid):
        for x, node in enumerate(row):
            #first compute neighbours from the particles binned in the last step
            #FIXME: perform the calculation necessary to get around the above fixme
            node.assignNeighbours()
            # print(f"Node at [{minX + Node.gridspacing*x},{minY + Node.gridspacing*y}] has {len(node.nList)}")

            #only consider computing defects if we are above a critical density in this node
            if (node.localDensity > Node.delta_density):
                currDefect = node.evalTopo()
                if currDefect != None: # if we didn't get the value 0 back, then we have a defect!
                    defect_container.append(currDefect)

    # return np.asarray([defect.getDefectData() for defect in defect_container])
    return defect_container

def singleFrame(input_filename : str,
                output_filename : str=None,
                particle_class : Type[Particle]=RodShapedBacterium):
    #construct infile name
    # input_filename = input_filename + str(step).zfill(5) + "." + outType
    step = int(re.findall("\d+",input_filename)[-1])
    print(f"step: {step}")

    # load data file into memory and store as a particle list
    print("Reading file " + input_filename + " for data.")
    try:
        data = pd.read_csv(input_filename,sep="\t")
        particle_list = intialiseElementsFromData(data,particle_class)
    except Exception as e:
        print("Failed to read data, terminating.")
        quit()

    totParticleNo = len(particle_list)
    print(f"Finished loading data. Total particle number: {totParticleNo}")

    # -------------------------- sort into layers ------------------------------
    maxZ=0
    # layerContainer = [] # lists of particles in each slice
    z_layers = round(maxZ / Node.layerheight) # number of layers we expect to use
    # for i in range(z_layers): # set up each layer's list
    #     layerContainer.append([])
    #print("Splitting into layers. Expecting " + str(z_layers) + " layers.")

    # No need to store all layers at once
    #loop through particle list and place them into appropriate list based on
    # for part in particle_list:
    #     scaledHeight = part.posZ / Node.layerheight
    #     targetLayer = int(scaledHeight) # the integer part only
    #
    #     if abs(scaledHeight-targetLayer) < EPSILON_LAYER: #if you're within epsilon of the target layer
    #         layerContainer[targetLayer].append(part) #move particle to the appropriate layer


    plane_loc_list = [0] # just to get started
    defect_container = [] # store any located defects

    for z_index in [z_layers]:

        plane_loc = z_index * Node.layerheight

        # Filter out all particles which interect the plane \pm eps
        # For now project onto xy plane
        layer = [particle for particle in particle_list
                 if cellPlaneIntersect(particle,plane_loc=plane_loc,
                 eps=Node.epsilon_layer,select=2)
                 ]

        print("There are {} particles in layer {} \pm {}".format(len(layer),
                                                                 plane_loc,
                                                                 Node.epsilon_layer)
                                                                 )
        defect_container.append(singleFrameFromCellData(particle_list=layer,
                                                        plane_loc=plane_loc)
                                                        )

    # flatten list
    defect_container = [defect for layer in defect_container for defect in layer]

    #------------------------- handle file writing -----------------------------

    #construct outfile path
    if output_filename is None:
        output_filename="./data/defectout_{:05d}.txt".format(step)

    open(output_filename, "w").close() #force clear the file
    with open(output_filename, "w") as file:
        #prepare file with headers
        file.write("id\tposX\tposY\tposZ\tcharge\torientation\n")
        #now loop through and dump defects
        if len(defect_container) > 0:
            for i, defect in enumerate(defect_container):
                defect.dump(i, file)
            print("\nDumped defects to " + output_filename + "\n")
        else:
            print(f"No defects found for step {step}!")

if __name__ == "__main__":
    #------------------- Handle Script Arguments -------------------------------
    # Example usage:
    # python3 ./FindDefects.py /path/to/frame/input_name.txt /path/to/frame/output_name.txt 1 10 0.5
    print( "Arguments:" )
    for arg in sys.argv:
        print( "\t" + arg )
    print()
    # path to data .txt file
    aDATANAME = sys.argv[1]

    # path to write dumped defects to
    aOUTNAME = sys.argv[2]

    # the layer height, acts as a characteristic scale for the script
    aLAYERHEIGHT = float(sys.argv[3])

    # How many "slices" to make in the x and y direction for the grid.
    # Should be comparible length scale to RodShapedBacterium length
    aGRIDSPACING = int(sys.argv[4])

    # how much to vary spacing when searching for defects
    aSEARCHSPACING = float (sys.argv[5])

    EPSILON_LAYER = 0.1
    EPSILON_LOOP = 0.05
    EPSILON_TOPOTHRESHOLD = 0.00001
    # R_NEIGHBOUR = 20
    # R_DENSITY = 20
    R_NEIGHBOUR = 1
    R_DENSITY = 1e-5
    # DELTA_DENSITY = 0.1
    DELTA_DENSITY = 1e-5
    N_LOOPMIN = 25

    # Set Node with these I/O and defect calculation hyperparameter data
    Node.setHyperParams(aLAYERHEIGHT,
                        aGRIDSPACING,
                        aSEARCHSPACING,
                        epsilon_layer=EPSILON_LAYER,
                        epsilon_loop=EPSILON_LOOP,
                        epsilon_topothreshold=EPSILON_TOPOTHRESHOLD,
                        r_neighbour=R_NEIGHBOUR,
                        r_density=R_DENSITY,
                        delta_density=DELTA_DENSITY,
                        n_loopmin=N_LOOPMIN
                        )

    # ---------------------- Search for Defects --------------------------------
    singleFrame(aDATANAME,aOUTNAME,particle_class=RodShapedBacterium)
