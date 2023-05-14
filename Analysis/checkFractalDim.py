# Standard libraries
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import axis,cm
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
fontprops = fm.FontProperties(size=10)
import numpy as np
import numpy.linalg as linalg
import scipy
import scipy.optimize as scopt
import pandas as pd
import re
import os
from pprint import pprint
from functools import reduce
from shapely.geometry import Polygon,Point,MultiPoint
from shapely.affinity import scale,rotate
from descartes import PolygonPatch
from matplotlib.collections import PatchCollection, LineCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from glob import glob
from functools import reduce
from scipy.spatial.ckdtree import cKDTree
from copy import copy

# third party modules
import imageio
from tqdm import tqdm
from skimage.morphology import skeletonize, medial_axis
from joblib import Parallel, delayed
import multiprocessing

# Custom modules
from DistributionFunctions import (getRadialDistributionFunction,
                                   getColonyDensity,
                                   computeColonyContour,
                                   getClusterContoursDescriptors,
                                   # boxCount,
                                   getContourCurvature,
                                   fitEllipse,
                                   checkBuckled,
                                   getAreaFraction,
                                   filterOutliers,
                                   getCurvatureGridMethod)
from FindDefects import intialiseElementsFromData
from visualise_infected_biofilm import VisInfectedBiofilm
from visualise_biofilm import setFigSize
from fastCellPlotting import addAllCellsToPlot
from phage_vis_class import Phage
from infected_vis_class import Infected
from RodShapedBacteria import RodShapedBacterium
from ChainingRodShapedBacteria import ChainingRodShapedBacterium
from SphericalBacteria import SphericalBacterium

# def boxCount(coords):
#     """
#     Find the fractal dimension by box counting
#     """
#     def makeSquare(mp,sl):
#         hsl=0.5*sl
#         a=np.array([-1,-1])*hsl
#         b=np.array([1,-1])*hsl
#         return Polygon([mp+a,mp+b,mp-a,mp-b,mp+a])
#
#     # diffs=np.diff(coords,axis=0)
#     # dists=np.linalg.norm(diffs,axis=1)
#     # mx_dist=np.max(dists) # The maximum length of a contour link
#
#     # Make the contour into a shapely object
#     # contourgon = Polygon([(v[0],v[1]) for v in coords])
#     points=MultiPoint([Point([v[0],v[1]]) for v in coords])
#     mx_dist=100
#     divisions = np.logspace(2,9,base=2,num=8,dtype=int)
#     num_boxes = []
#     side_lengths = []
#     left, down, right, up = points.bounds
#     for nn in divisions:
#         # start=time.perf_counter()
#
#         origin = np.array([left,down])
#         bounding_side_length=np.max([right-left,up-down])
#         xs,box_side = np.linspace(0,bounding_side_length,num=nn,
#                                   endpoint=True,retstep=True
#                                   )
#         xx,yy = np.meshgrid(xs,xs)
#         cx=xx.flatten()
#         cy=yy.flatten()
#
#         midpoints=np.array([mp for mp in zip(cx,cy)])
#         box_tree=cKDTree(midpoints+origin[np.newaxis,:])
#         point_tree=cKDTree(coords)
#         # neighbour_boxes=point_tree.query_ball_point(midpoints,
#         #                                             r=0.5*np.sqrt(2)*box_side
#         #                                             )
#         neighbour_pnts=box_tree.query_ball_tree(point_tree,
#                                                 r=0.5*np.sqrt(2)*box_side
#                                                 )
#         boxes=[]
#         mps=[]
#         for ii in range(len(neighbour_pnts)):
#             mp=midpoints[ii]+origin
#             box=makeSquare(mp,box_side)
#             for jj in neighbour_pnts[ii]:
#                 p=Point(coords[jj])
#                 if box.contains(p):
#                     boxes.append(box)
#                     mps.append(mp)
#                     break
#         mp_check_tree=cKDTree(mps)
#         res=mp_check_tree.sparse_distance_matrix(mp_check_tree,
#                                                  max_distance=np.sqrt(2)*box_side,
#                                                  output_type='coo_matrix')
#         res.setdiag(1e10)
#         if np.min(res.data)==0:
#             print("Error! Some boxes were too far away from each other")
#             quit()
#         side_lengths.append(box_side)
#         num_boxes.append(len(boxes))
#
#         # end=time.perf_counter()
#         # print(f"{nn} divisions ({box_side=}) took {end-start} seconds")
#         if nn==16:
#             example_boxes=copy(boxes)
#
#         # fig,ax = plt.subplots(1,1,figsize=setFigSize(247))
#         # ax.scatter(coords[:,0],coords[:,1],lw=0.7,zorder=3)
#         # for box in boxes:
#         #     ax.plot(*box.exterior.xy,'k',alpha=0.5,lw=0.7,zorder=1)
#         # ax.axis("scaled")
#         # plt.show()
#         # plt.close()
#
#     # The gradient of the linear fit is the fractal dimension
#     linear_fit = lambda x,m,c: m*x+c
#     ln_num = np.log(num_boxes)
#     ln_sides = -np.log(side_lengths)
#     popt,pcov = scopt.curve_fit(linear_fit,
#                                 ln_sides,
#                                 ln_num
#                                 )
#     # plt.plot(ln_sides,ln_num,
#     #          marker='o',mfc='w',mec='k'
#     #          )
#     # plt.plot(ln_sides,linear_fit(ln_sides,*popt),'k--')
#     # plt.show()
#     # print(popt)
#     return popt[0],example_boxes

fname='../Analysis/koch_snowflake.png'

# Create an ImageJ2 gateway with the newest available version of ImageJ2.
import imagej
ij = imagej.init()
# Load an image.
image_url = 'https://www.researchgate.net/profile/Jack-Tuszynski/publication/311503809/figure/fig1/AS:437008762970112@1481202416370/An-image-of-the-Koch-snowflake-a-fractal-with-fractal-dimension-d-126-From.png'
jimage = ij.io().open(image_url)

# Convert the image from ImageJ2 to xarray, a package that adds
# labeled datasets to numpy (http://xarray.pydata.org/en/stable/).
image = ij.py.from_java(jimage)

# Display the image (backed by matplotlib).
ij.py.show(image, cmap='gray')

quit()

img=imageio.imread(fname)
select=np.where(img[:,:,1]==0)
nx,ny,_=img.shape
xs,ys=np.mgrid[0:nx,0:ny]
X=xs[select]
Y=ys[select]
coords=np.array([cc for cc in zip(X,Y)])

# imgplot = plt.imshow(img)
# koch_snowflake=Polygon(coords)
# fig,ax=plt.subplots()
# ax.axis("scaled")
# plt.plot(*koch_snowflake.exterior.coords.xy)
# plt.scatter(coords[:,0],coords[:,1],c='r')
# ax.set_xlim(X.min(),X.max())#
# ax.set_ylim(Y.min(),Y.max())
# plt.show()
# quit()
fd,boxes=boxCount(coords)
print(f"{fd=}")
