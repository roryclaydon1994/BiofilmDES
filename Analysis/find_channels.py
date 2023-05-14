"""

    Find channels from BiofilmPhageDES projection outputs

    This script contiains the channel location algorithm. It will assume the
    following files are on the path:
    RodShapedBacteria.py
    ChainingRodShapedBacteria.py
    visualise_biofilm.py

    Example usage is provided at the bottom under direct invocation

    Expected usage:

    python3 find_channels.py

Author: Rory Claydon

"""

# Standard modules
import pandas as pd
import matplotlib as mpl
import re
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.pyplot import axis, cm
from matplotlib.colors import ListedColormap
from scipy.interpolate import griddata
import numpy as np
import glob
import os.path
import sys
import multiprocessing as mp
import pandas
from numba import jit, cuda
from tqdm import tqdm
from scipy.spatial.transform import Rotation
from joblib import Parallel, delayed
import multiprocessing as mp
from pprint import pprint

# Connected-component analysis
import scipy.sparse as sps
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.ndimage import label
import networkx as nx
from scipy.spatial.ckdtree import cKDTree
from scipy.spatial import ConvexHull,convex_hull_plot_2d
from skimage import measure
from skimage.draw import ellipsoid

# third party
try:
    from polylabel import polylabel
except Exception as e:
    print(e)

# User defined classes
from RodShapedBacteria import RodShapedBacterium
from ChainingRodShapedBacteria import ChainingRodShapedBacterium

def pointToLineSeg(p,r1,r2):
    """
        Find the shortest distance from a point to a line segment in d dims

        Based on Dan Sunday's implementation
        http://geomalgorithms.com/a02-_lines.html

        Parameters:
            p: float, array like
                the location vector of the point
            r1: float, array like
                the start of the line segment
            r2: float, array like
                the end of the line segment

        Returns:
            dist: float
                shortest distance to the line segment
    """
    v = r2-r1  # direction vector
    w = p-r1   # point to line start

    vw = v.dot(w) # projection onto line segment, negative implies obtuse angle
    vv = v.dot(v) # squared length of line seg

    if vw < 0:
        # point is before the line start
        return np.linalg.norm(w)
    elif vv < vw:
        # point is beyond the end of the line
        return np.linalg.norm(p-r2)
    else:
        # the right angle from the point to the line is within the segment
        mu = vw / vv  # fraction along seg where the right angle to p intersects
        r = r1 + mu * v
        return np.linalg.norm(p-r)

def isPointInCell(p,cell,thresh=None,dim="2d"):
    """
        Find if a vector lies within a cell

        Parameters:
            p: float, array like
                the vector of the point to check
            cell: instance of cell class
                the cell to check if the point is inside
            thresh: float
                threshold distance to be considered to be inside the cell.
                If None, then use cell radius
        Returns:
            inside: bool
                True if point is in this cell
    """

    if thresh == None:
        thresh=0.5*cell.diameter

    if dim=="2d":
        cta = cell.tangent_axes
    else:
        cta=[0,1,2]
    cell_upper_pole = cell.rcm[cta] + 0.5*cell.length*cell.ori[cta]
    cell_lower_pole = cell.rcm[cta] - 0.5*cell.length*cell.ori[cta]
    return pointToLineSeg(p,cell_lower_pole,cell_upper_pole) < thresh

def findPointsInCloud(a,b,threshold,box_widths,dim="2d"):
    ab = b-a
    cm = 0.5*(a+b)
    L = 2*threshold+np.linalg.norm(ab)
    R = 2*threshold+1
    try:
        dl = 0.9*min(box_widths) # take min as rotate to actual orientation later
    except Exception as e:
        dl = 0.9*box_widths # take min as rotate to actual orientation later

    x = np.arange(-0.5*L,0.5*L+dl,dl)
    y = np.arange(-0.5*R,0.5*R+dl,dl)
    if dim=="2d":
        X,Y = np.meshgrid(x,y)
        coords=np.asarray([X.flatten(),Y.flatten()])
        phi = np.arctan2(ab[1],ab[0])
        c,s = np.cos(phi),np.sin(phi)
        rot=np.array([[c,-s],[s,c]])
        return (cm[:,np.newaxis]+rot.dot(coords)).T
    elif dim=="3d":
        X,Y,Z = np.meshgrid(x,y,y)
        coords=np.asarray([X.flatten(),Y.flatten(),Z.flatten()])
        r = generateRotation(ab)
        t_pnts=(cm[np.newaxis,:]+r.apply(coords.T))
        return t_pnts
    else:
        print(f"{dim} is invalid dimension")
        quit()

def getIndices(point,box_widths,Ls):
    indices=np.floor((0.5*Ls+point)/box_widths)
    return [int(ii) for ii in indices]

def findBiofilmLimits(grid):
    """
    Find the parts of grid which mark the top of the biofilm.

    If grid[i,j,k]=1, this is a vacancy.

    If grid[i,j,k]=0, this is occupied.
    """

    grid_pts = grid.shape

    s = np.ones((3,3))  # Allow all connections
    for kk in range(grid_pts[0]):
        print("---k---")
        print(kk)
        l,n=label(grid[kk],s)   # get lables and connected segments
        print("---n---")
        print(n)
        print("---g---")
        outer_space_indices = np.where(l<=1)
        grid[kk][outer_space_indices] = 0
        print(grid[kk])
        print("---l---")
        print(l)
        print("-----")

def getCavities3D(cell_list,ax_lims,box_widths,threshold=1,mode="DBSCAN",
                  use_checkpoint=None):
    """
        Locate channel structures in 3D.

        Parameters:
            cell_list: list of susceptible instances
                cells to determine cavities between
            ax_lims: 3,2 matrix
                The limits in which to search for cavities
            box_widths: array [3], float
                side lengths of the grid box (x,y,z)
            threshold: float
                how far outside a cell is still considered "inside" cell.
                Essentially smoothes out the cell
            mode: str
                "DBSCAN" or "CCA". Method of locating connected components

        Returns:
            cavity_data: data_frame
                the columns are "x","y","z","k" for the coordinates in 3 space
                and the cluster identifier
    """

    """
    TODO: move the fill in space with bacteria component to separate function
    TODO: try voxel representation of cavities
    """
    ax_rng = np.diff(ax_lims,axis=1).flatten()
    grid_pts=1+np.floor(ax_rng/box_widths)
    Ls=grid_pts*box_widths
    axes  = [np.arange(-0.5*(Ls[ii]-box_widths[ii]),0.5*Ls[ii],box_widths[ii])
                for ii in range(3)]

    grid_coords = np.meshgrid(*axes,sparse=False,indexing='ij')

    if use_checkpoint == True:
        grid = np.loadtxt("test/grid.csv")
        cavity_data=pd.read_csv("test/cavity_data.csv",sep='\t')
        grid=grid.reshape([int(ii) for ii in grid_pts[::-1]])
        # print(grid.shape)
        print("loaded!")
    else:

        grid=np.ones([int(pp) for pp in grid_pts[::-1]])

        for cell in tqdm(cell_list):
            cta = cell.tangent_axes
            a = cell.rcm + 0.5*(cell.length+1)*cell.ori
            b = cell.rcm - 0.5*(cell.length+1)*cell.ori
            trial_points=findPointsInCloud(a,b,threshold,box_widths,"3d")
            # try this with numba
            # plt.close()
            # cell.addElementToPlot(ax,ax_rng=ax_rng)
            # fig,ax=plt.subplots(1,1,subplot_kw={'projection': "3d"})
            missed_points=0
            for trial_point in trial_points:
                # Generate a list of (x,y) points to check for occupancy
                if isPointInCell(trial_point,cell,threshold,dim="3d"):
                    # If a grid cell centre is within threshold consider it occupied
                    # print(trial_point)
                    indices=getIndices(trial_point,box_widths,Ls=ax_rng)
                    # ax.scatter(*trial_point,s=10,fc='r')
                    # print(f"setting {indices} to 0")
                    try:
                        grid[indices[2],indices[1],indices[0]]=0
                    except Exception as e:
                        missed_points+=1
                        # print(e)
                        # print(f"{trial_point} is indexed {indices} for {0.5*Ls} in {grid.shape[::-1]}")
                        # print(f"skipping {trial_point}")
            # print(f"skipped {missed_points}/{len(trial_points)} at {cell.rcm}")
                # else:
                #     ax.scatter(*trial_point,s=10,fc='k')


            # ax.plot(*[[a[ii],b[ii]] for ii in range(3)])
            # plt.show()
            # quit()
        print("saving test/grid.csv and test/cavity_data.csv")
        np.savetxt("test/grid.csv",grid.flatten())

        # Remove anything 'outside' biofilm
        # findBiofilmLimits(grid)

        cs = np.argwhere(grid > 0)
        cavity_data=detectBlobs(cs,eps=1,mode=mode)
        cavity_data.to_csv("test/cavity_data.csv",sep='\t')

    cavities_num=len(cavity_data['k'].unique())
    # print(f"Found {cavities_num} clusters")
    largest_cavities=cavity_data['k'].value_counts().nlargest(5)
    # print(largest_cavities)

    return cavity_data, grid, grid_coords

    # print("------------")
    # print("Finding convex hull of cells")
    # cell_grid = 1 - grid
    # points=np.argwhere(cell_grid>0)
    # cell_hull=ConvexHull(points)
    # print("------------")
    # cs = np.argwhere(grid > 0)
    # print(cs.shape)
    #
    # print("------------")
    # print("Clustering cavities")
    # cavity_data=detectBlobs(cs,eps=1,mode=mode)

    # fig = plt.figure(figsize=(10,10))
    # ax = fig.add_subplot(111, projection='3d')
    # for ii in range(0,grid.shape[0],5):
    #
    #     # if grid[ii].min() != 0:
    #     #     continue
    #
    #     # cd = cavity_data.loc[cavity_data['x']==ii].to_numpy()
    #
    #     vertices=cell_hull.vertices[points[cell_hull.vertices][:,0]==ii]
    #     origin=np.array([min(axes[0]),min(axes[1]),min(axes[2])])
    #     print(f"origin at {origin}")
    #     z=origin[2]+box_widths[2]*ii
    #     print(f"At z {z}")
    #
    #     # hull_coords=origin[np.newaxis,:]+box_widths[np.newaxis,:2]*points[vertices][:,[2,1]]
    #     # hull2d = ConvexHull(hull_coords)
    #     # uncomment to plot the convex hull
    #     # convex_hull_plot_2d(hull2d,ax=ax)
    #
    #     print("Finding 2d contours")
    #     contours = measure.find_contours(grid[ii])
    #     n=len(contours)
    #     viridis = cm.get_cmap('viridis', n)
    #     newcolors = viridis(np.linspace(0, 1, n))
    #     # blank = np.array([0, 0, 0, 0])
    #     # newcolors[0:2, :] = blank
    #     newcmp = ListedColormap(newcolors)
    #     for contour,colour in zip(contours,np.random.permutation(newcolors)):
    #         if len(contour)>10:
    #             coords = origin[np.newaxis,:2]+box_widths[np.newaxis,:2]*contour[:,[1,0]]
    #             # com = np.average(coords,axis=0)
    #             vis_centre = polylabel([coords])
    #             # ax.scatter(vis_centre[0],vis_centre[1],z,s=20,c=[colour],ec='k')
    #             ax.plot(coords[:,0],coords[:,1],z,"-",lw=2,color=colour)
    #
    #     # ax.plot(hull_coords[hull2d.vertices,0],hull_coords[hull2d.vertices,1 ], 'r--', lw=2)
    #     ax.set_xlim([-0.5*ax_rng[0],0.5*ax_rng[0]])
    #     ax.set_ylim([-0.5*ax_rng[1],0.5*ax_rng[1]])
    #     # ax.axis("scaled")
    #
    #     # cavities_num=len(cavity_data.loc[cavity_data['x']==ii]['k'].unique())
    #     # print(f"Found {cavities_num} clusters")
    # plt.show()
    # quit()



    # for ii in range(grid.shape[0]):
    #     # if grid[ii].min() != 0:
    #     #     continue
    #     if grid[ii].min() != 0:
    #         continue
    #
    #     cd = cavity_data.loc[cavity_data['x']==ii].to_numpy()
    #     fig = plt.figure(figsize=(10,10))
    #     ax = fig.add_subplot(111, projection=None)
    #
    #     # n=2
    #     # viridis = cm.get_cmap('viridis', n)
    #     # newcolors = viridis(np.linspace(0, 1, n))
    #     # blank = np.array([0, 0, 0, 0])
    #     # newcolors[0, :] = blank
    #     # newcmp = ListedColormap(newcolors)
    #     # ax.imshow(grid[ii],origin='lower',
    #     #           extent=[-0.5*ax_rng[0],0.5*ax_rng[0],
    #     #                   -0.5*ax_rng[1],0.5*ax_rng[1]],
    #     #           cmap=newcmp
    #     #           )
    #
    #     for cell in cell_list:
    #         cell.addElementToPlot(ax,ax_rng=ax_rng[0])
    #
    #     vertices=cell_hull.vertices[points[cell_hull.vertices][:,0]==ii]
    #     origin=np.array([min(axes[0]),min(axes[1])])
    #     print(f"origin at {origin}")
    #     hull_coords=origin[np.newaxis,:]+box_widths[np.newaxis,:2]*points[vertices][:,[2,1]]
    #     hull2d = ConvexHull(hull_coords)
    #     # uncomment to plot the convex hull
    #     # convex_hull_plot_2d(hull2d,ax=ax)
    #
    #     print("Finding 2d contours")
    #     contours = measure.find_contours(grid[ii])
    #     n=len(contours)
    #     viridis = cm.get_cmap('viridis', n)
    #     newcolors = viridis(np.linspace(0, 1, n))
    #     # blank = np.array([0, 0, 0, 0])
    #     # newcolors[0:2, :] = blank
    #     newcmp = ListedColormap(newcolors)
    #     for contour,colour in zip(contours,np.random.permutation(newcolors)):
    #         if len(contour)>10:
    #             coords = origin[np.newaxis,:]+box_widths[np.newaxis,:2]*contour[:,[1,0]]
    #             # com = np.average(coords,axis=0)
    #             vis_centre = polylabel([coords])
    #             ax.scatter(vis_centre[0],vis_centre[1],s=20,c=[colour],ec='k')
    #             ax.plot(coords[:,0],coords[:,1],"-",lw=2,color=colour)
    #
    #     # ax.plot(hull_coords[hull2d.vertices,0],hull_coords[hull2d.vertices,1 ], 'r--', lw=2)
    #     ax.set_xlim([-0.5*ax_rng[0],0.5*ax_rng[0]])
    #     ax.set_ylim([-0.5*ax_rng[1],0.5*ax_rng[1]])
    #     ax.axis("scaled")
    #     plt.show()
    #     plt.close()
    #
    #     cavities_num=len(cavity_data.loc[cavity_data['x']==ii]['k'].unique())
    #     print(f"Found {cavities_num} clusters")
    #     # quit()
    #
    # return cavity_data, grid_coords

def getCavities(cell_list,ax_rng,box_width=0.1,threshold=1):
    # 2D for now, and in xy plane to check it works
    grid, XX, YY = get2DOccupations(cell_list,ax_rng,box_width,threshold)

    s = np.ones((3,3))  # Allow all connections
    l,n=label(grid,s)   # get lables and connected segments

    print("---n---")
    print(n)
    print("---g---")
    print(grid)
    print("---l---")
    print(l)
    print("-----")
    return l,n,XX,YY

def get2DOccupations(cell_list,ax_rng,box_width=0.1,threshold=1):
    """
        Find a coarse grid of whether points are occupied by cells
    """
    # 2D for now, and in xy plane to check it works
    assert(cell_list[0].projection=="xy")
    N=int(1+ax_rng//box_width)
    L=N*box_width
    x=np.arange(-0.5*(L-box_width),0.5*L,box_width)
    XX,YY=np.meshgrid(x,x,sparse=False,indexing='ij')

    print("----------")
    print(f"box_width: {box_width} threshold: {threshold} N: {N}")

    grid=np.ones((N,N))    # 1 for no occupation, 0 else

    for cell in tqdm(cell_list):
        cta = cell.tangent_axes
        a = cell.rcm[cta] + 0.5*(cell.length+1)*cell.ori[cta]
        b = cell.rcm[cta] - 0.5*(cell.length+1)*cell.ori[cta]
        trial_points=findPointsInCloud(a,b,threshold,box_width)
        for trial_point in trial_points:
            # Generate a list of (x,y) points to check for occupancy
            if isPointInCell(trial_point,cell,threshold):
                # If a grid cell centre is within threshold consider it occupied
                indices=getIndices(trial_point,box_width,Ls=L)
                try:
                    grid[indices[0],indices[1]]=0
                except Exception as e:
                    pass

    return grid, XX, YY
