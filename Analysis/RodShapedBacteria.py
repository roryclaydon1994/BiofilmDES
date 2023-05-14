# Standard libraries
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import math
import numpy.linalg as linalg
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.patches import Wedge
from shapely.geometry import Polygon,Point
from shapely.ops import unary_union

# third party
from tqdm import tqdm
# import geopandas as gpd
import pandas as pd
from scipy.spatial.ckdtree import cKDTree
from scipy.ndimage.filters import minimum_filter,maximum_filter
import networkx as nx

# from sklearn.cluster import DBSCAN
from numba import jit,njit,prange,set_num_threads

# If changing to conda we can use this
# import cudf
# import cuspatial

def getCharge(sorted_loop):
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

        perp_to_d2 = np.array([dir2[1],-dir2[0]])

        det = np.dot(dir1,perp_to_d2)
        ang = math.atan2(abs(det), d1_dot_d2)

        if ang > 0.5*np.pi: #flip as necessary
            dir2 = -dir2
            det = -det
            d1_dot_d2 = -d1_dot_d2

        # add the angle to the topo counter
        phi += np.sign(det)*math.atan2(abs(det),d1_dot_d2)

    return 0.5*phi/np.pi # return the topo counter

class RodShapedBacterium(object):
    """
        Control visualisation of susceptible in a biofilm

        Author: Rory Claydon

    """
    elastic_moduli=1
    sus_vis_radius_factor=1
    depletion_radius=None
    colour_fn=lambda th: cm.get_cmap('hsv')(getNematicAngle(th)/(np.pi))
    smoothCutOff=6 # When smoothing cells, this is the max dist from rcm to consider
    cell_hash = {}

    def __init__(self,cell_id,length,radius,
                      pos_x,pos_y,pos_z,
                      ori_x,ori_y,ori_z,
                      virtual_center_x=0.0,
                      virtual_center_y=0.0,
                      virtual_center_z=0.0,
                      force_x=0.0,
                      force_y=0.0,
                      force_z=0.0,
                      torque_x=0.0,
                      torque_y=0.0,
                      torque_z=0.0,
                      neighbours=None,
                      projection="xy"):

        self.cell_id          = int(cell_id)    # cell id
        self.length           = length          # pole to pole length
        self.radius           = radius          # spherocylinder radius
        self.diameter         = 2*radius        # spherocylinder diameter
        self.virtual_center_x = virtual_center_x # virtual centre x coordinate for a given interaction
        self.virtual_center_y = virtual_center_y
        self.virtual_center_z = virtual_center_z
        self.force_x          = force_x          # force x coordinate on cell
        self.force_y          = force_y
        self.force_z          = force_z
        self.torque_x         = torque_x         # torque x coordinate on cell
        self.torque_y         = torque_y
        self.torque_z         = torque_z
        self.colour           = (30/255,130/255,76/255,1) # default element colour

        try:
            self.neighbours = [ int(id) for id in neighbours.split(',')
                                               if id != '' ]
        except Exception as e:
            # print(e)
            try:
                self.neighbours = id
            except Exception as e:
                # print(e)
                print("Error loading neighbour list")
                self.neighbours = None
                quit()

        self.rcm = np.array([pos_x,pos_y,pos_z])

        self.ori  = np.array([ori_x,ori_y,ori_z])

        self.virt = np.array([virtual_center_x,virtual_center_y,virtual_center_z])

        self.frce = np.array([force_x,force_y,force_z])

        self.tau  = np.array([torque_x,torque_y,torque_z])

        self.theta = np.arctan2(ori_y,ori_x)
        self.nematic_angle = getNematicAngle(self.theta)

        self.nematic_ori=self.getNematicOrientation()

        self.energy=0

        # Given a projection determine the in and out of plane axes
        # This creates attricutes tangent_axes and normal_axes
        # See documentation updateProjection for details
        # This will be run every time a particle is initialised which is
        # inefficient. TO DO: Move to calling functions
        self.updateProjection(projection)
        projection = projection

    def __str__(self):
        cell_str = "cell {:d}, L:{:3.3f} d:{:3.3f} at ({:3.3f},{:3.3f},{:3.3f})"
        return cell_str.format(self.cell_id,self.length,self.diameter,*self.rcm)

    @classmethod
    def makeHashList(cls,cells):
        # Convert to hash list for fast lookup
        cls.cell_hash={ cell.cell_id : cell for cell in cells }

    @classmethod
    def updateProjection(cls,new_projection):
        """
            Project the cell onto a new plane. Works by defining the class
            variables tangent_axes and normal_axes.

            Parameters:
                new_projection: str
                    The plane onto which to project
                cls: class
                    class to update the class variables of
            Effect:
                updates self.tangent_axes:
                    The plane tangent axes i.e. [0,1] for axis 0 and 1 (x and y)
                    intended used length[tangent_axes] gives length vector
                    in this plane
                updates self.normal_axes
                    The axis the normal points in i.e. length[normal_axes]
                    is the out of plane length
        """
        cls.projection = new_projection
        tangent_axes, normal_axes = setProjection(new_projection)
        cls.tangent_axes = tangent_axes
        cls.normal_axes = normal_axes

    def getArea(self):
        """
            Compute the area of the projection of this particle onto the specified
            plane
        """
        in_plane_length = self.ori[self.tangent_axes]*self.length
        in_plane_length = np.linalg.norm(in_plane_length)
        radius = 0.5*self.diameter
        return in_plane_length*self.diameter + np.pi*radius**2

    def getVolume(self):
        """
            Compute the 3D volume of this particle
        """
        radius = 0.5*self.diameter
        return np.pi*self.length*radius**2 + (4*np.pi*radius**3)/3

    def getInPlaneDirector(self):
        """
            Find the normalised in plane director
        """
        director2D = self.ori[self.tangent_axes]
        return getUnitVec(director2D)

    def getPolygon(self,radius=None):
        """
            Get a shapely representation of the cell
        """
        if radius==None:
            radius=self.radius

        # center of mass + half full length * tangent to bacteria for top pole
        center_top    = self.rcm + 0.5*self.length*self.ori
        center_bottom = self.rcm - 0.5*self.length*self.ori

        self.updateProjection('xy')

        # Perp vec for drawing the polygon in the plane for the body of the rod
        perp_vector = np.array([self.ori[self.tangent_axes[1]],
                                -self.ori[self.tangent_axes[0]]
                                ],
                                dtype = float)
        if linalg.norm(perp_vector)>0:
            perp_vector /= linalg.norm(perp_vector)

        # Make the main body to fill in the colour
        points = [center_top[self.tangent_axes]   - radius*perp_vector,
                  center_top[self.tangent_axes]   + radius*perp_vector,
                  center_bottom[self.tangent_axes]+ radius*perp_vector,
                  center_bottom[self.tangent_axes]- radius*perp_vector]

        main_body = Polygon(points)
        head = Point(center_top).buffer(radius)
        tail = Point(center_bottom).buffer(radius)
        rep = unary_union([head,main_body,tail])
        return rep

    def addElementToPlot(self,ax,colour=None,projection="xy",ax_rng=20,show_id=False,zorder=1,ec='k'):
        """
            Draw a susceptible on an axis
            Parameters:
                ax: axis handle
                    The axis to draw the particle on
                colour: colour
                    Specify fill colour for this element
                projection: str
                    determines onto which plane to project the 3D biofilm
                    options: "xy","xz","yz"
                annotate: str
                    annotate element with this string
                zorder: int
                    layer level in plot
            Returns:
                perp_vector: float array, 3 elements
                    The vector perpendicular to the orientaion
        """

        if colour==None:
            colour=self.colour

        plot_radius=0.5*self.diameter*RodShapedBacterium.sus_vis_radius_factor

        cell=self.getPolygon(plot_radius)

        ax.fill(*cell.exterior.xy,fc=colour,ec=ec,alpha=0.8,linewidth=5/ax_rng,zorder=zorder)

        # # center of mass + half full length * tangent to bacteria for top pole
        # center_top    = self.rcm + 0.5*self.length*self.ori
        # center_bottom = self.rcm - 0.5*self.length*self.ori
        #
        # self.updateProjection(projection)
        #
        # # Place the circles at the poles for the ends of the rod
        # for center in [center_top,center_bottom]:
        #     circle = plt.Circle(center[self.tangent_axes],
        #                         radius=plot_radius,
        #                         fc=colour,
        #                         edgecolor='k',
        #                         linewidth=5/ax_rng)
        #     ax.add_patch(circle)
        #
        # # Perp vec for drawing the polygon in the plane for the body of the rod
        # perp_vector = np.array([self.ori[self.tangent_axes[1]],
        #                         -self.ori[self.tangent_axes[0]]
        #                         ],
        #                         dtype = float)
        # if linalg.norm(perp_vector)>0:
        #     perp_vector /= linalg.norm(perp_vector)
        #
        # # Make the main body to fill in the colour
        # points = [center_top[self.tangent_axes]   - plot_radius*perp_vector,
        #           center_top[self.tangent_axes]   + plot_radius*perp_vector,
        #           center_bottom[self.tangent_axes]+ plot_radius*perp_vector,
        #           center_bottom[self.tangent_axes]- plot_radius*perp_vector]
        # line = plt.Polygon(points, closed=True, facecolor=colour,
        #                            edgecolor=None,linewidth=5/ax_rng)
        # ax.add_patch(line)
        #
        # # Make the edges to colour to give the cell a nice boundary
        # points = [center_top[self.tangent_axes]   - plot_radius*perp_vector,
        #           center_bottom[self.tangent_axes]- plot_radius*perp_vector]
        # line = plt.Polygon(points, closed=False, facecolor=colour,
        #                            edgecolor="k",linewidth=5/ax_rng)
        # ax.add_patch(line)
        #
        # # Make the edges to colour to give the cell a nice boundary
        # points = [center_top[self.tangent_axes]   + plot_radius*perp_vector,
        #           center_bottom[self.tangent_axes]+ plot_radius*perp_vector]
        # line = plt.Polygon(points, closed=False, facecolor=colour,
        #                            edgecolor="k",linewidth=5/ax_rng)
        # ax.add_patch(line)
        #
        # # Testing depletion interaction
        # # if RodShapedBacterium.depletion_radius!=None:
        # #     Rd=RodShapedBacterium.depletion_radius+0.5*self.diameter
        # #     theta1=self.theta-np.pi/2
        # #     theta2=theta1+np.pi
        # #     for center in [center_top,center_bottom]:
        # #         arc = Wedge(center[self.tangent_axes],
        # #                     r=Rd,
        # #                     theta1=theta1*180/np.pi,
        # #                     theta2=theta2*180/np.pi,
        # #                     fc='none',edgecolor='k',
        # #                     ls='--',linewidth=5/ax_rng)
        # #         ax.add_patch(arc)
        # #         theta1+=np.pi
        # #         theta2=theta1+np.pi
        # #
        # #     # Make the edges to colour to give the cell a nice boundary
        # #     points = [center_top[self.tangent_axes]   - Rd*perp_vector,
        # #               center_bottom[self.tangent_axes]- Rd*perp_vector]
        # #     line = plt.Polygon(points, closed=False, facecolor=None,
        # #                                edgecolor="k",
        # #                                ls='--',linewidth=5/ax_rng)
        # #     ax.add_patch(line)
        # #
        # #     # Make the edges to colour to give the cell a nice boundary
        # #     points = [center_top[self.tangent_axes]   + Rd*perp_vector,
        # #               center_bottom[self.tangent_axes]+ Rd*perp_vector]
        # #     line = plt.Polygon(points, closed=False, facecolor=None,
        # #                                edgecolor="k",
        # #                                ls='--',linewidth=5/ax_rng)
        # #     ax.add_patch(line)
        #
        if show_id==True:
            pnt=cell.representative_point()
            ax.text(pnt.x,pnt.y,f"{self.cell_id}")
        # return perp_vector

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):

        patches=[]

        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        radius=0.5*height
        diameter=2*radius
        tangent_axes=[0,1]
        colour=self.colour
        lw=0.1

        # center of mass + half full length * tangent to bacteria for top pole
        ori=np.array([1.0,0.0])
        rcm=np.array([x0+0.5*width,y0+0.5*height])
        center_top    = rcm + 0.5*width*ori
        center_bottom = rcm - 0.5*width*ori

        # Place the circles at the poles for the ends of the rod
        for center in [center_top,center_bottom]:
            circle = mpatches.Circle(center,
                                     radius=radius,
                                     fc=colour,edgecolor='k',lw=lw,
                                     transform=handlebox.get_transform())
            patches.append(circle)
            #
            # arc = Arc(center,height=diameter,width=diameter,angle=0,theta1=0,theta2=90,lw=lw)
            # patches.append(arc)

        # Perp vec for drawing the polygon in the plane for the body of the rod
        perp_vector = np.array([ori[tangent_axes[1]],-ori[tangent_axes[0]]])
        # print(perp_vector)
        # perp_vector /= linalg.norm(perp_vector)
        # print(perp_vector)

        # Make the main body to fill in the colour
        points = [center_top[tangent_axes]   - radius*perp_vector,
                  center_top[tangent_axes]   + radius*perp_vector,
                  center_bottom[tangent_axes]+ radius*perp_vector,
                  center_bottom[tangent_axes]- radius*perp_vector]
        line = plt.Polygon(points, closed=True, facecolor=colour,
                                   edgecolor=None,lw=lw,
                                   transform=handlebox.get_transform())
        patches.append(line)

        # Make the edges to colour to give the cell a nice boundary
        points = [center_top[tangent_axes]   - radius*perp_vector,
                  center_bottom[tangent_axes]- radius*perp_vector]
        line = plt.Polygon(points, closed=False, facecolor=colour,
                                   edgecolor="k",lw=lw,
                                   transform=handlebox.get_transform())
        patches.append(line)

        # Make the edges to colour to give the cell a nice boundary
        points = [center_top[tangent_axes]   + radius*perp_vector,
                  center_bottom[tangent_axes]+ radius*perp_vector]
        line = plt.Polygon(points, closed=False, facecolor=colour,
                                   edgecolor="k",lw=lw,
                                   transform=handlebox.get_transform())
        patches.append(line)

        patch_collection = PatchCollection(patches,match_original=True)
        handlebox.add_artist(patch_collection)
        return patch_collection

    def addHertzianInteractionToPlot(self,ax):
        """
            Draw a virtual spheres, forces and torques on an axis
            Parameters:
                ax: axis handle
                    The axis containing the elements on which to add the forces
        """
        circle = plt.Circle((self.virt[0],self.virt[1]),
                             radius=0.5*self.diameter, facecolor="None", edgecolor='r')
        ax.add_patch(circle)
        ax.scatter(self.virt[0],self.virt[1],40,facecolor="r",zorder=10)
        # ax.text(self.virt[0],self.virt[1],"vc=({},{},{})".format(*self.virt))

        force_norm = self.frce / np.linalg.norm(self.frce)
        ax.arrow(x=self.virt[0],y=self.virt[1],
                 dx=force_norm[0],dy=force_norm[1],
                 width=0.05,facecolor='k')

        # Make sure axis will show everything
        ax.scatter(1.1*(self.virt[0] + force_norm[0]),
                   1.1*(self.virt[1] + force_norm[1]),0.0001)

        ax.text(1.1*(self.virt[0] + force_norm[0]),
                1.1*(self.virt[1] + force_norm[1])-0.15,
                "f=({:3.3f},{:3.3f},{:3.3f})".format(*self.frce))

    def getSmoothCellRho(self,r,sigma=1):
        """
            Get a smoothed density profile for this cell

            Parameters:
                r: array, float
                    position at which the profile is observed
                sigma: float
                    smoothing length
            Returns:
                rho: float
                    density at this location
        """
        delta_r = r[self.tangent_axes] - self.rcm[self.tangent_axes]
        if np.linalg.norm(delta_r)>RodShapedBacterium.smoothCutOff:
            return 0

        tang = self.getInPlaneDirector()
        perp = np.array([tang[1],-tang[0]])
        perp/=np.linalg.norm(perp)

        eff_hlf_L = 0.5*(self.length+self.diameter) # effective half length
        eff_hlf_W = 0.5*self.diameter               # effective half width

        smooth_main = 0.5*( np.tanh( (tang.dot(delta_r)+eff_hlf_L) / sigma )
                           -np.tanh( (tang.dot(delta_r)-eff_hlf_L) / sigma )
                           )
        smooth_perp = 0.5*( np.tanh( (perp.dot(delta_r)+eff_hlf_W) / sigma )
                           -np.tanh( (perp.dot(delta_r)-eff_hlf_W) / sigma )
                           )
        return smooth_main*smooth_perp

    def getSmoothCellQ(self,r):
        """
            Find this cell's contribution to the Q tensor at position r
        """
        tang = self.getInPlaneDirector()
        return self.getSmoothCellRho(r)*( 2*np.outer(tang,tang)-np.identity(2) )

    def getSmoothCellEnergy(self,r):
        """
            Find this cell's contribution to the energy at position r
        """
        return self.getSmoothCellRho(r)*self.energy

    def getEnergy(self):
        """
            The only energy for the base case is the Hertzian energy
        """
        self.hertzian_energy=0
        for cell_b in self.contacts:
            assert( np.isclose(self.radius,cell_b.radius) )
            h=(self.radius+cell_b.radius)-getCellDist(self,cell_b)
            h_energy = (2/5) * RodShapedBacterium.elastic_moduli * h**(2/5)
            self.hertzian_energy+=h_energy
        self.hertzian_energy*=0.5 # account for over counting
        self.energy+=self.hertzian_energy

    def getNematicOrientation(self):
        return np.array([np.cos(self.nematic_angle),np.sin(self.nematic_angle)])

    @staticmethod
    def getContactEnergies(cells):
        if hasattr(cells[0],'contacts')==False:
            RodShapedBacterium.getNeighbourLists(cells)

        for cell in cells:
            cell.getEnergy()

        return [cell.hertzian_energy for cell in cells]

    @staticmethod
    def computeSmoothedRho(cell_list,r):
        """
            Find a global density profile by smoothing out the cells
        """
        rho_list = [ cell.getSmoothCellRho(r) for cell in cell_list ]
        return sum(rho_list)

    @staticmethod
    def computeSmoothedEnergy(cell_list,r):
        """
            Find a global density profile by smoothing out the cells
        """
        energy_list = [ cell.getSmoothCellEnergy(r) for cell in cell_list ]
        return sum(energy_list)

    @staticmethod
    def getCellGridNeighbours(grid_coords,cell_list):
        """
            Bin the grid points to the nearest cell locations
        """

        grid_tree = cKDTree(grid_coords)

        X = np.array([ cell.rcm[:2] for cell in cell_list ])

        # result = grid_tree.query_ball_point(
        #     X,
        #     r=RodShapedBacterium.smoothCutOff,
        #     workers=-1
        #     )
        cell_tree = cKDTree(X)
        result = cell_tree.query_ball_tree(
            grid_tree,
            r=RodShapedBacterium.smoothCutOff
            )

        max_neighbours=max([ len(r) for r in result ])
        nat_result=-np.ones((len(result),max_neighbours),dtype=int)
        for ii,nlist in enumerate(result):
            for jj,nn in enumerate(nlist):
                nat_result[ii,jj]=nlist[jj]
        cell_data = RodShapedBacterium.toNative(cell_list)

        return nat_result, cell_data

    @staticmethod
    def computeDensity(cell_list,dr=0.5,fname='rho.npy'):
        """
            Find the smoothed density profile on a grid
        """
        #----------------------- Create grid ---------------------------
        grid = RodShapedBacterium.createGrid(cell_list,dr)
        xv,yv = grid

        try:
            density_grid=np.load(fname)
        except Exception as e:
            print(e)
            grid_coords = np.array([xv.flatten(),yv.flatten()]).T

            #------------------ Find cell neighbours on grid ----------------------
            nat_result,cell_data = RodShapedBacterium.getCellGridNeighbours(grid_coords,cell_list)

            #--------------------- Find density on grid --------------------
            # Column and row element numbers
            ny,nx = xv.shape
            density_grid = np.zeros(ny*nx)

            print("Launch threads")
            try:
                set_num_threads(12)
            except Exception as e:
                set_num_threads(1)

            findDensityGrid(nat_result,cell_data,density_grid,grid_coords,
                            RodShapedBacterium.smoothCutOff
                            )

            density_grid=density_grid.reshape((ny,nx))
            np.save(fname,density_grid)

        return density_grid, grid

    @staticmethod
    def computeEnergyGrid(cell_list,dr=0.5):
        """
            Find the smoothed density profile on a grid
        """

        #--------------------------- Create grid -------------------------------
        xv,yv = RodShapedBacterium.createGrid(cell_list,dr)

        # Column and row element numbers
        ny,nx = xv.shape
        energy_grid = np.zeros(ny*nx)

        #------------------------- Find energy on grid ----------------------------

        grid_coords = np.array([xv.flatten(),yv.flatten()]).T
        cell_tree = cKDTree(grid_coords)

        X = np.array([ cell.rcm[:2] for cell in cell_list ])

        result = cell_tree.query_ball_point(
            X,
            r=RodShapedBacterium.smoothCutOff,
            workers=1
            )
        max_neighbours=max([ len(r) for r in result ])
        nat_result=-np.ones((len(result),max_neighbours),dtype=int)
        for ii,nlist in enumerate(result):
            for jj,nn in enumerate(nlist):
                nat_result[ii,jj]=nlist[jj]
        cell_data = RodShapedBacterium.toNative(cell_list)
        cell_energy = np.array([ cell.energy for cell in cell_list ])

        findEnergyGrid(nat_result,cell_data,cell_energy,energy_grid,grid_coords,
                       RodShapedBacterium.smoothCutOff
                       )

        energy_grid=energy_grid.reshape((ny,nx))
        return energy_grid, [xv,yv]

    @staticmethod
    def computeSmoothedQ(cell_list,r):
        """
            Find the Q-tensor using the smoothed cell density
        """
        Qs = [ cell.getSmoothCellQ(r) for cell in cell_list ]
        Q = np.sum(Qs,axis=0)
        return Q

    @staticmethod
    def computeLocalCharge(cell_list,r,dr=0.01):
        """

        """
        r_x_dr = np.array([dr,0])
        Q_x_p = RodShapedBacterium.computeSmoothedQ(cell_list,r+r_x_dr)
        Q_x_m = RodShapedBacterium.computeSmoothedQ(cell_list,r-r_x_dr)
        Q_dx = ( Q_x_p-Q_x_m ) / ( 2*dr )

        r_y_dr = np.array([0,dr])
        Q_y_p = RodShapedBacterium.computeSmoothedQ(cell_list,r+r_y_dr)
        Q_y_m = RodShapedBacterium.computeSmoothedQ(cell_list,r-r_y_dr)
        Q_dy = ( Q_y_p-Q_y_m ) / ( 2*dr )

        Qxx_dx = Q_dx[0,0]
        Qxy_dx = Q_dx[0,1]
        Qxx_dy = Q_dy[0,0]
        Qxy_dy = Q_dy[0,1]

        q = ( Qxx_dx*Qxy_dy - Qxy_dx*Qxx_dy ) / (2*np.pi)
        return q

    @staticmethod
    def computeChargeDensity(cell_list,dr=0.1,fname="Q_grid.npy"):
        """
            Find the local charge density from the Q-tensor
        """

        #--------------------------- Create grid -------------------------------
        grid = RodShapedBacterium.createGrid(cell_list,dr)
        xv,yv = grid
        ny,nx = xv.shape

        try:
            Q=np.load(fname)
        except Exception as e:
            print(e)
            grid_coords = np.array([xv.flatten(),yv.flatten()]).T

            #------------------ Find cell neighbours on grid ----------------------
            nat_result,cell_data = RodShapedBacterium.getCellGridNeighbours(grid_coords,cell_list)

            #------------------------- Find Q on grid ------------------------------
            Q_grid = np.zeros((ny*nx,2,2))

            print("Launch threads")
            set_num_threads(12)
            findQGrid(nat_result,cell_data,Q_grid,grid_coords,
                      RodShapedBacterium.smoothCutOff
                      )

            Q=Q_grid.reshape((ny,nx,2,2))
            np.save(fname,Q)

        #--------------- Try to compute q with central derivatives -------------
        # Only calculate on the inner grid for now
        q = np.zeros((ny,nx))
        print("finding q")
        Q_dx = np.gradient(Q,axis=0,edge_order=2)
        Q_dy = np.gradient(Q,axis=1,edge_order=2)
        q = (Q_dx[:,:,0,0]*Q_dy[:,:,0,1]-Q_dx[:,:,0,1]*Q_dy[:,:,0,0])/(2*np.pi)

        return q, Q, grid

    @staticmethod
    def computeChargeFromLoop(r,cells):
        """
            Find the local charge density from the Q-tensor
            Based on Tim's work
        """
        def getIntersectionMidpoint(cc,loop):
            line=cc.getPolygon().intersection(loop)
            if hasattr(line,'geoms'):
                line=max([ll for ll in line.geoms],key=lambda x: x.length)
            p=line.interpolate(0.5,normalized=True)
            return p.coords[0]

        def getTheta(cc,loop):
            p=getIntersectionMidpoint(cc,loop)
            return  math.atan2( *(p-r)[[1,0]] )

        eps=0.05*RodShapedBacterium.smoothCutOff
        success=False
        while eps<RodShapedBacterium.smoothCutOff:
            circle=Point(r).buffer(eps)
            loop=circle.exterior
            # in_circle=[ cell.getPolygon().intersection(circle) for cell in cells ]
            selected_cells=[ cell for cell in cells
                             if cell.getPolygon().intersects(loop)
                             ]
            num_loop=len(selected_cells)
            # dens=np.sum([poly.area for poly in in_circle])/circle.area
            charge=0
            # if num_loop>3 and dens>0.25:
            if num_loop>5:
                gt=lambda cc: getTheta(cc,loop)
                sorted_loop = sorted(selected_cells,
                                     key=gt
                                     )
                charge=getCharge(sorted_loop)
                if abs(abs(charge)-0.5)<0.1:
                    success=True
                    # fig,ax=plt.subplots(1,1)
                    # for cell in selected_cells:
                    #     cell.colour="r"
                    #     ax.scatter(*getIntersectionMidpoint(cell,loop),zorder=4)
                    # for cell in cells:
                    #     cell.addElementToPlot(ax)
                    # for poly in in_circle:
                    #     ax.fill(*poly.exterior.xy,fc='b',ec='k',alpha=0.8,zorder=2)
                    # ax.plot(*loop.xy)
                    # ax.axis("scaled")
                    # ax.set_title(f"{charge=:3.3e}")
                    # fig.savefig(f"../GeneratedOutput/tmp/{r=}_{charge=:3.3f}.pdf",
                    #     format='pdf',
                    #     bbox_inches='tight',
                    #     transparent=True
                    #     )
                    break
            eps+=eps
        # print(f"{r} {charge} {success}")
        return charge,success

    @staticmethod
    def clusterDefects(X):
        tree=cKDTree(X)
        result=tree.query_ball_point(X,r=5)
        reduced_points=[]
        checked_points=[]

        for ii in range(len(X)):
            points_to_check=result[ii]
            cluster=[]
            for point in result[ii]:
                if point not in checked_points:
                    checked_points.append(point)
                    cluster.append(X[point])
            if cluster:
                reduced_points.append(np.mean(cluster,axis=0))
        return np.array(reduced_points)

    @staticmethod
    def filterDefects(q_field,grid):
        """
            Filter defects
        """
        X,Y=grid
        filters=[minimum_filter,maximum_filter]
        defect_positions=[]
        for filter in filters:
            res=(q_field==minimum_filter(q_field,3,mode='constant',cval=0.0))
            locs=np.where(res==1)
            qs=q_field[locs]
            xs=X[locs]
            ys=Y[locs]
            thresh=np.abs(qs)>0.05*np.max(np.abs(qs))
            dfcts=RodShapedBacterium.clusterDefects(
                np.array([xs[thresh],ys[thresh]]).T
                )
            defect_positions.append(dfcts)
        return defect_positions

    @staticmethod
    def findLocalDefects(X,cells):
        """
            Given a set of potential defect points and a list of cells, try to find
            if the defect in the vicinity is a real defect classified by if the
            strength is greater than 0.2
        """
        cell_data = RodShapedBacterium.toNative(cells)
        cell_tree=cKDTree(cell_data[:,:2])
        point_tree=cKDTree(X)

        # Find the cell neighbour indices of the suggested defects
        nn = point_tree.query_ball_tree(
            cell_tree,
            r=RodShapedBacterium.smoothCutOff
            )

        dxs=[]
        charges=[]
        for ii,x in enumerate(X):
            ncs=[cells[jj] for jj in nn[ii]]
            q,success=RodShapedBacterium.computeChargeFromLoop(x,ncs)
            if success:
                charges.append(q)
                dxs.append(x)
        return np.array(charges),np.array(dxs)

    @staticmethod
    def createGrid(cell_list,dr=0.1):
        """
            Make a grid given the cell list and resolution of points
        """
        RodShapedBacterium.updateProjection('xy')
        tang_axes = RodShapedBacterium.tangent_axes

        # Find grid boundaries
        cell_rcms = [ cell.rcm[tang_axes] for cell in cell_list ]
        cell_rcms = np.asarray(cell_rcms)

        # Add on increase by the maximum half length to include ends of cells
        min_x,min_y=np.min(cell_rcms,axis=0) - 3.5
        max_x,max_y=np.max(cell_rcms,axis=0) + 3.5

        # Ensure all cells are enclosed
        x = np.arange(min_x,max_x+dr,dr)
        y = np.arange(min_y,max_y+dr,dr)

        xv,yv = np.meshgrid(x,y,sparse=False,indexing='xy')

        return xv,yv

    @staticmethod
    def getCOM(cell_list):
        com=np.array([0.0,0.0,0.0])
        mass=0
        for cell in cell_list:
            volume=cell.getVolume()
            com+=volume*cell.rcm
            mass+=volume
        return com/mass

    @staticmethod
    def getOverlaps(cells):
        RodShapedBacterium.getNeighbourLists(cells)
        overlaps=[]
        for cc in cells:
            for nn in cc.contacts:
                if nn.cell_id<cc.cell_id:
                    toverlap=(cc.radius+nn.radius)-getCellDist(cc,nn)
                    overlaps.append(toverlap)
        return overlaps

    @staticmethod
    def checkContact(cell_a,cell_b):
        return getCellDist(cell_a,cell_b)<(cell_a.radius+cell_b.radius)

    @staticmethod
    def getNeighbourLists(cells,r=6,find_contacts=True):
        N=len(cells)
        X = np.array([ cell.rcm[:2] for cell in cells ])
        cell_tree = cKDTree(X)
        result = cell_tree.query_ball_point(X,r=r)

        if find_contacts==False:
            for ii in range(N):
                cells[ii].neighbours=[]
                for jj in result[ii]:
                    if jj < ii:
                        cells[ii].neighbours.append(cells[jj])
                        cells[jj].neighbours.append(cells[ii])
        else:
            for ii in range(N):
                cells[ii].neighbours=[]
                for jj in result[ii]:
                    if jj < ii:
                        max_sep=(
                            cells[ii].radius+cells[jj].radius+
                            0.5*(cells[ii].length+cells[jj].length)
                            )
                        if np.linalg.norm(cells[ii].rcm-cells[jj].rcm)<max_sep:
                            cells[ii].neighbours.append(cells[jj])
                            cells[jj].neighbours.append(cells[ii])

            for cc in cells:
                cc.contacts=[]
                for nn in cc.neighbours:
                    if RodShapedBacterium.checkContact(cc,nn)==True:
                        cc.contacts.append(nn)

    def getCellQ(self):
        """
            This could be faster as only need the first row
        """
        return 2*np.outer(self.ori[:-1],self.ori[:-1])-np.identity(2)

    @staticmethod
    def findMicrodomains(cells,colour_fn,threshold_angle=3):
        RodShapedBacterium.getNeighbourLists(cells) # Populate neighbour lists if in contact
        print("Try connected component analysis")

        # G=nx.Graph([(1,2),(2,3),(4,5)])
        # nx.draw(G)
        # plt.show()
        # for s in nx.connected_components(G):
        #     print(s)

        def isPaired(cc,nn):
            dp=np.dot(cc.ori,nn.ori)
            denom = np.linalg.norm(cc.ori)*np.linalg.norm(nn.ori)
            x=dp/denom
            if abs(x)>=1:
                x=1
            theta=np.arccos( x )
            theta=min(theta,np.pi-theta)
            return theta<=threshold_angle*np.pi/180

        pairs=[(cc.cell_id,nn.cell_id) for cc in cells
               for nn in cc.contacts if isPaired(cc,nn)
               ]
        # N=len(cells)
        G=nx.Graph(pairs)
        # G.add_nodes_from([ cell.cell_id for cell in cells ])
        largest_cc = max(nx.connected_components(G), key=len)
        print(f"{len(largest_cc)=}")

        clusters={}
        # cell_hash={ cell.cell_id : cell for cell in cells }
        for ii,s in enumerate(nx.connected_components(G)):
            cid=ii+1
            clusters[cid]=[ RodShapedBacterium.cell_hash[id] for id in s  ]

        for key,cells in clusters.items():
            cc=cells[0]
            for nn in cells:
                nn.colour=colour_fn(cc.theta)

        # for cc in tqdm(cells):
        #     G.add_node(cc.cell_id)
        #     for nn in cc.contacts:
        #         mat.append([cc.cell_id,nn.cell_id,1])
        return clusters
        # quit()

        clusters={}
        current_cluster=0
        used_ids=[]

        for cc in tqdm(cells):
            if cc.cell_id not in used_ids:
                current_cluster+=1
                used_ids.append(cc.cell_id)
                try:
                    clusters[f"{current_cluster}"].append(cc)
                except Exception as e:
                    clusters[f"{current_cluster}"]=[cc]
                cells_to_check=cc.contacts
                cc.colour=colour_fn(cc.theta)
                while cells_to_check:
                    nn=cells_to_check[0]
                    available_cell=(
                        isPaired(cc,nn) and
                        nn.cell_id not in used_ids
                        )
                    if available_cell:
                        # print(cc.cell_id,nn.cell_id)
                        nn.colour=colour_fn(cc.theta)
                        clusters[f"{current_cluster}"].append(nn)
                        used_ids.append(nn.cell_id)
                        cells_to_check+=nn.contacts
                    else:
                        cells_to_check.remove(nn)
        print(f"Found {current_cluster} clusters")
        # for key,cl_cells in clusters.items():
        #     print(f"cluster {key}")
        #     cls=""
        #     for cc in cl_cells:
        #         cls+=f"{cc.cell_id} "
        #     print(cls)
        return clusters

    def getNativeData(self):
        return [*self.rcm[:2], *self.ori[:2], self.length, self.diameter]

    @staticmethod
    def toNative(cells):
        return np.array([ cell.getNativeData() for cell in cells ])
"""
    Define utilities
"""
def setProjection(new_projection):
    """
        Get the plane normal and tangent axes given a projection

        Parameters:
            new_projection: str
                The plane onto which to project
        Returns:
            tangent_axes: array like, int
                The plane tangent axes i.e. [0,1] for axis 0 and 1 (x and y)
                intended used length[tangent_axes] gives length vector
                in this plane
            normal_axes: int
                The axis the normal points in i.e. length[normal_axes] is
                the out of plane length
    """
    if new_projection == "xy":
        tangent_axes = [0,1]
        normal_axes = 2
    elif new_projection == "xz":
        tangent_axes = [0,2]
        normal_axes = 1
    elif new_projection == "yz":
        tangent_axes = [1,2]
        normal_axes = 0
    return tangent_axes,normal_axes

def getUnitVec(vec):
    """
        Find the unit vector given an input vector
        Parameters:
            vec: array-like, float
                vector to find the unit vector of
        Return:
            unit_vec: array-like, float
                unit vector
    """
    return vec / np.linalg.norm(vec)

def getLineToLineParams(x_1,a_1,x_2,a_2):
    """
    Move this to a geometry module later
    """
    def clampRoot(root):
        if (root>=1):
            return 1
        elif (root<=-1):
            return -1
        else:
            return root

    x_21 = x_1 - x_2 # vector from x2 to x1
    a =  a_1.dot(a_1)
    b = -a_1.dot(a_2)
    c =  a_2.dot(a_2)
    d =  a_1.dot(x_21)
    e = -a_2.dot(x_21)
    delta = a*c - b*b   # \equiv modulus( a_1 \cross a_2 )

    # === Handle degenerate cases ===
    if ( a==0 and c==0 ):
        s_star=0
        t_star=0
    elif ( a==0 and c>0 ):
        s_star=0
        t_star=clampRoot(-e/c)
    elif ( c==0 and a>0 ):
        t_star=0
        s_star=clampRoot(-d/a)
    else:
        # === Handle line segments ===
        # Handle parallel line segments
        if (delta==0):
           s_star = 0.5*( clampRoot( (c-e)/b )+clampRoot( -(c+e)/b ) )
           t_star = clampRoot(-(e+b*s_star)/c)
        else:
            t_star = clampRoot( ( b*d - a*e ) / delta )
            s_star_tmp = (-b*t_star - d) / a
            s_star = clampRoot( s_star_tmp )
            if ( (s_star_tmp < -1) or (s_star_tmp > 1) ):
                t_star = clampRoot( ( -b*s_star - e ) / c )
    return s_star,t_star

def getCellDist(cell_a,cell_b):
    a_1 = 0.5*cell_a.ori*cell_a.length
    a_2 = 0.5*cell_b.ori*cell_b.length
    x_1 = cell_a.rcm
    x_2 = cell_b.rcm
    s,t=getLineToLineParams(x_1=x_1,a_1=a_1,x_2=x_2,a_2=a_2)
    c1=x_1+s*a_1
    c2=x_2+t*a_2
    dist=np.linalg.norm(c1-c2)
    return dist

def getNematicAngle(alpha):
    if np.abs(alpha)>np.pi:
        print(f"alpha is too big! {alpha:2.3f}")
        quit()

    if alpha>=0:
        return alpha

    if alpha<0:
        return np.pi+alpha

@njit
def getSmoothCellRhoNat(cell_data,r,smoothCutOff):
    """
        Get a smoothed density profile for this cell using native types

        This will only work in 2D

        Parameters:
            cell_data:
                stuff
            r: array, float
                position at which the profile is observed
            sigma: float
                smoothing length
        Returns:
            rho: float
                density at this location
    """
    sigma=1
    rcm=cell_data[0:2]
    tang=cell_data[2:4]
    length=cell_data[4]
    diameter=cell_data[5]

    delta_r = r - rcm
    if np.linalg.norm(delta_r)>smoothCutOff:
        return 0

    perp = np.array([tang[1],-tang[0]])
    perp/=np.linalg.norm(perp)

    eff_hlf_L = 0.5*(length+diameter) # effective half length
    eff_hlf_W = 0.5*diameter          # effective half width

    smooth_main = 0.5*( np.tanh( (tang.dot(delta_r)+eff_hlf_L) / sigma )
                       -np.tanh( (tang.dot(delta_r)-eff_hlf_L) / sigma )
                       )
    smooth_perp = 0.5*( np.tanh( (perp.dot(delta_r)+eff_hlf_W) / sigma )
                       -np.tanh( (perp.dot(delta_r)-eff_hlf_W) / sigma )
                       )
    return smooth_main*smooth_perp

@njit
def getSmoothCellEnergyNat(cell_data,cell_energy,r,smoothCutOff):
    return getSmoothCellRhoNat(cell_data,r,smoothCutOff)*cell_energy

@njit
def getSmoothCellQNat(cell_data,r,smoothCutOff):
    """
        Find this cell's contribution to the Q tensor at position r
    """
    tang=cell_data[2:4]
    density_cont = getSmoothCellRhoNat(cell_data,r,smoothCutOff)
    return density_cont*( 2*np.outer(tang,tang)-np.identity(2) )

@njit(parallel=True)
def findDensityGrid(result,cell_data,density_grid,grid_coords,
                   smoothCutOff):
    # for cc,pnt in enumerate(result):
    n,m = result.shape
    for ii in prange(n):
        for jj in range(m):
            if result[ii,jj]>=0:
                pnt=result[ii,jj]
                density_grid[pnt] += getSmoothCellRhoNat(cell_data[ii],
                                                         grid_coords[pnt],
                                                         smoothCutOff
                                                         )

@njit(parallel=True)
def findQGrid(result,cell_data,Q_grid,grid_coords,
                   smoothCutOff):
    # for cc,pnt in enumerate(result):
    n,m = result.shape
    for ii in prange(n):
        for jj in range(m):
            if result[ii,jj]>=0:
                pnt=result[ii,jj]
                Q_grid[pnt,:,:] += getSmoothCellQNat(cell_data[ii],
                                                     grid_coords[pnt],
                                                     smoothCutOff
                                                     )

@njit
def findEnergyGrid(result,cell_data,cell_energy,energy_grid,grid_coords,
                   smoothCutOff):
    # for cc,pnt in enumerate(result):
    n,m = result.shape
    for ii in range(n):
        for jj in range(m):
            if result[ii,jj]>=0:
                pnt=result[ii,jj]
                energy_grid[pnt] += getSmoothCellEnergyNat(cell_data[ii],
                                                           cell_energy[ii],
                                                           grid_coords[pnt],
                                                           smoothCutOff
                                                           )

@njit
def test(result,cell_data,cell_energy,energy_grid,grid_coords,smoothCutOff):
    return smoothCutOff+cell_data[3]
