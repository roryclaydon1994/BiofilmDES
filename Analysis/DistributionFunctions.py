"""
    Analysis script which provides distribution function analysis to help
    understand colony morphology.

"""

# Standard modules
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import axis, cm
from matplotlib import ticker
from matplotlib.colors import ListedColormap
from skimage import measure
import scipy
import scipy.optimize as scopt
from scipy.interpolate import UnivariateSpline, splev, splprep
from copy import copy
import time
from matplotlib.collections import PatchCollection, LineCollection
from descartes import PolygonPatch
import scipy.stats as stats
import os

#Third party modules
import shapely.speedups
shapely.speedups.enable()
import shapely.vectorized
from shapely.affinity import scale,rotate,translate
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union
from scipy.spatial.ckdtree import cKDTree
from scipy.ndimage.filters import maximum_filter
from skimage.morphology import skeletonize, medial_axis
from skimage.measure import EllipseModel
from numba import jit,njit,prange,set_num_threads
from tqdm import tqdm

# Custom modules
import fastCellPlotting as fcp

# params={'font.family'           :'serif',
#         'text.usetex'           : True,
#         'axes.titlesize'        : 10,
#         'axes.labelsize'        : 10,
#         'xtick.labelsize'       : 10,
#         'ytick.labelsize'       : 10,
#         'legend.frameon'        : False,
#         'lines.linewidth'       : 1.5,
#         'legend.handlelength'   : 1.5,
#         'legend.columnspacing'  : 0.75,
#         'legend.handletextpad'  : 0.5,
#         # 'legend.title_fontsize' : None,
#         'legend.fontsize'       : 10,
#         'font.size'             : 10,
#         'legend.columnspacing'  : 1.0,
#         "axes.formatter.useoffset":False,
#         "text.latex.preamble"   : [r"\usepackage[detect-all,locale=DE]{siunitx}",r"\usepackage[T1]{fontenc}",r"\usepackage{amsmath}"]
#        }
# mpl.rcParams.update(params)

def makeAnnulus(c,r1,r2):
    big=Point(c).buffer(r2)
    small=Point(c).buffer(r1)
    return big.difference(small)

try:
    from polylabel import polylabel
except Exception as e:
    print(e)

# Custom modules
from RodShapedBacteria import RodShapedBacterium
from AG43RodShapedBacteria import AG43RodShapedBacterium
from ChainingRodShapedBacteria import ChainingRodShapedBacterium
from visualise_biofilm import setFigSize
from find_channels import get2DOccupations
from visualise_infected_biofilm import VisInfectedBiofilm

def filterOutliers(arr):
    """
        Given an array, filter the outliers using the 1.5 IQR rule
    """
    if len(arr)>1:
        iqr=stats.iqr(arr)
        q1=np.percentile(arr,q=25)
        q3=np.percentile(arr,q=75)
        idx=(q1-1.5*iqr<arr) & (arr<q3+1.5*iqr)
        accepted=arr[ idx ]
        outliers=arr[ ~idx ]

        print(
            f"{iqr=}",
            f"{q1=}",
            f"{q3=}",
            f"{idx=}",
            f"{accepted=}",
            f"{outliers=}"
            )
    else:
        accepted=arr
        outliers=np.empty(1)
    return accepted,outliers

def getRadialDistributionFunction(cells):
    """
        Compute radial distribution function

        Parameters:
            cells: list of cell instances

        Returns:
            rs: float, array
                centres of the bins
            probs: float, array
                probability of being between r[i] and r[i+1]

    """
    com=RodShapedBacterium.getCOM(cells)
    print(f"{com=}")
    radial_dist=np.zeros((len(cells),1))
    for ii,cell in enumerate(cells):
        if np.isclose(cell.rcm[2],0)==False:
            print("We expect 2d colonies for this colony")
            quit()
        radial_dist[ii]=np.linalg.norm(cell.rcm-com)
    counts, bins = np.histogram(radial_dist,bins=50,density=True)
    dr=bins[1]-bins[0]
    rs=bins[:-1]+0.5*dr
    probs=counts/(2*np.pi*rs*dr)
    return rs,probs

def computeColonyContourHoles(cells):
    """
        Given a list of cells, find the largest contour polygon in terms of area.
    """
    cellgons=[ cell.getPolygon(radius=cell.radius) for cell in cells ]
    try:
        springs = [ spring.buffer(1e-3*cell.radius) for cell in cells
                    for spring in cell.getSpringGons(radiusfactor=0.7*cell.radius)
                    ]
        cellgons.extend(springs)
    except Exception as e:
        # print(e)
        pass
    cell_repr = unary_union(cellgons)

    try:
        out_contour=max(cell_repr.geoms,key=lambda x: x.area)
    except Exception as e:
        out_contour=cell_repr
    return out_contour

def computeColonyContour(cells,add_links=True,ret_gon=False):
    """
        Find the colony outer colony (if one exists)

        Parameters:
            cells: list of Cell instances
            smoothing: float
                multiply the radius by this

        Returns:
            coords: (N,2) numpy array
                x,y coordinates of the outer colony
    """
    # ax_lim = np.max(VisInfectedBiofilm.getAxisLimits(cells)) + 5
    # grid,_,_ = get2DOccupations(cell_list=cells,ax_rng=ax_lim,
    #                             box_width=box_width,threshold=threshold)
    # contours = measure.find_contours(grid)
    # longest_contour = max(contours,key=len)
    # coords = -0.5*ax_lim+box_width*longest_contour
    sm_fac=1
    cellgons=[ cell.getPolygon(radius=sm_fac*cell.radius) for cell in cells ]
    if add_links:
        try:
            # if not ChainingRodShapedBacterium.cell_hash:
            #     print("here")
                # ChainingRodShapedBacterium.makeHashList(cells)
            # springs=ChainingRodShapedBacterium.getAllSpringGons(cells)
            local_hash={ cell.cell_id : cell for cell in cells }
            springs = [ cell.getHeadSpringGon(local_hash) for cell in cells ]
            springs = [ spring for spring in springs if spring!=None ]
            # uncomment for springs for ag43
            # springs = [ spring.buffer(1e-3*cell.radius) for cell in cells
            #             for spring in cell.getSpringGons(radiusfactor=0.7*cell.radius)
            #             ]
            cellgons.extend(springs)
        except Exception as e:
            print(e)
            quit()
    cell_repr = unary_union(cellgons)

    # plt.close('all')
    # fig, ax = plt.subplots(1,1)
    max_length=0
    try:
        try:
            for geom in cell_repr.geoms:
                # ax.plot(*geom.exterior.xy)
                contour=geom
                if contour.length>max_length:
                    out_contour=contour
                    max_length=contour.length
        except Exception as e:
            # print(e)
            out_contour=cell_repr
    except Exception as e:
        print(e)
        quit()

    # print(f"{cell_repr.geom_type=}")

    # ax.plot(*out_contour.exterior.xy,'k--',lw=0.5)
    # for cc in cells:
    #     cc.addElementToPlot(ax)
    # ax.set_aspect('equal', 'datalim')
    contour_coords=np.array([(z[0],z[1]) for z in zip(*out_contour.exterior.xy)])
    # for cell in cells:
    #     if Polygon(contour_coords).contains(cell.getPolygon(radius=0.95*cell.radius))==False:
    #         print(f"error! {cell.cell_id} not in contour at {cell.rcm}")
    #         plt.close('all')
    #         RodShapedBacterium.getNeighbourLists(cells)
    #         fig,ax=plt.subplots(1,1)
    #         ax.plot(*out_contour.exterior.xy,'k--',lw=0.5)
    #         RodShapedBacterium.sus_vis_radius_factor=0.7
    #         cell.addElementToPlot(ax,ax_rng=2)
    #         for nc in cell.contacts:
    #             nc.addElementToPlot(ax)
    #         gonext=cell.getPolygon(radius=sm_fac*cell.radius).exterior
    #         l,d,r,u=gonext.bounds
    #         ax.plot(*gonext.xy)
    #         # try:
    #         #     springs = [ spring.buffer(1e-3*cell.radius) for cell in cells
    #         #                 for spring in cell.getSpringGons(radiusfactor=0.7*cell.radius)
    #         #                 ]
    #         #     patches = [ PolygonPatch(polygon=gon) for gon in springs ]
    #         #     ax.add_collection(PatchCollection(patches,match_original=True,zorder=1))
    #         # except Exception as e:
    #         #     print(e)
    #         ax.axis("scaled")
    #         xls=ax.set_xlim(l-cell.radius,r+cell.radius)
    #         yls=ax.set_ylim(d-cell.radius,u+cell.radius)
    #         print(f"Ax lims: {xls}, {yls}")
    #         fig.savefig(f"../GeneratedOutput/example_outside_contour_{cell.cell_id=}.pdf",
    #                     format="pdf",bbox_inches='tight',
    #                     transparent=True)
    #         plt.show()
    #         quit()
    if ret_gon:
        return out_contour
    else:
        return contour_coords

def getColonyDensity(cells,coords):
    """
        Find the density of the colony contained within the outer contour

        Parameters:
            cells: list of Cell instances
            coords: (N,2) float numpy array or None
                coords in x,y columns of the colony contour. If None, calculate
                this

        Returns:
            density: float
                density of cells contained in outer contour
            contour: shapely polygon
                extracted colony outer contour
    """
    cell_area=0
    for cell in cells:
        cell_area+=cell.getArea()
    try:
        contourgon = Polygon([(v[0],v[1]) for v in coords])
    except Exception as e:
        print(e)
        coords = computeColonyContour(cells)
        contourgon = Polygon([(v[0],v[1]) for v in coords])

    density=cell_area/contourgon.area
    return density,contourgon

def fitEllipse(coords):
    ell = EllipseModel()
    ell.estimate(coords)
    return ell.params

def getColonyEllipse(coords):
    params=fitEllipse(coords)
    xc, yc, a, b, theta = params
    circ=Point([xc,yc]).buffer(1)
    ellipse=rotate(scale(circ,0.5*a,0.5*b),theta,use_radians=True)
    return ellipse

def getClusterContoursDescriptors(clusters,cluster_contours,colony_contour,
                                  fig_name=None):
    """
    Find the average microdomain size
    """
# find the centre of the colony and the rough radius

    ellipse=getColonyEllipse(colony_contour)

    # plt.close()
    # plt.plot(*ellipse.exterior.xy)

    areas=np.zeros(len(clusters))
    interior_areas=[]
    for ii,(key,coords) in enumerate(cluster_contours.items()):
        p = Polygon([(v[0],v[1]) for v in coords])
        areas[ii]=p.area
        if p.intersects(ellipse):
            # plt.plot(*p.exterior.xy,color='r')
            interior_areas.append(p.area)
    #     else:
    #         plt.plot(*p.exterior.xy,color='k',linestyle=':')
    # plt.show()
    print(f"mean interior: {np.mean(interior_areas)} pm {np.std(interior_areas)}")

    if fig_name!=None:
        fig,ax=plt.subplots()
        ax.plot(*ellipse.exterior.xy,lw=1.5,color='b',ls='--')
        for ii,(key,coords) in enumerate(cluster_contours.items()):
            p = Polygon([(v[0],v[1]) for v in coords])
            if p.intersects(ellipse):
                plt.plot(*p.exterior.xy,color='r',lw=1,alpha=0.8)
                # ax.add_patch(PolygonPatch(p,fc='None',ec='r',alpha=0.8))
            else:
                plt.plot(*p.exterior.xy,color='k',lw=1,alpha=0.4)
        ax.axis('off')
        ax.axis('scaled')
        fig.savefig(fig_name,bbox_inches='tight',format='pdf',transparent=True)
        # plt.show()

        # plt.hist(areas,density=True,alpha=0.7,label='all')
        # plt.hist(interior_areas,density=True,alpha=0.7,label='interior')
        # plt.axvline(np.mean(interior_areas),color='m',alpha=0.7,label='mean interior')
        # plt.axvline(np.percentile(areas,q=75),color='k',alpha=0.7,label='q3 all')
        # plt.legend()
        # # plt.show()
        # plt.close('all')
    return np.mean(interior_areas)

def boxCount(coords):
    """
    Find the fractal dimension by box counting
    """
    def makeSquare(mp,sl):
        hsl=0.5*sl
        a=np.array([-1,-1])*hsl
        b=np.array([1,-1])*hsl
        return Polygon([mp+a,mp+b,mp-a,mp-b,mp+a])

    diffs=np.diff(coords,axis=0)
    dists=np.linalg.norm(diffs,axis=1)
    mx_dist=np.max(dists) # The maximum length of a contour link

    # Make the contour into a shapely object
    contourgon = Polygon([(v[0],v[1]) for v in coords])
    divisions = np.logspace(2,9,base=2,num=8,dtype=int)
    num_boxes = []
    side_lengths = []
    for nn in divisions:
        # start=time.perf_counter()
        left, down, right, up = contourgon.bounds
        origin = np.array([left,down])
        bounding_side_length=np.max([right-left,up-down])
        xs,box_side = np.linspace(0,bounding_side_length,num=nn,
                                  endpoint=True,retstep=True
                                  )
        xx,yy = np.meshgrid(xs,xs)
        cx=xx.flatten()
        cy=yy.flatten()

        midpoints=np.array([mp for mp in zip(cx,cy)])
        box_tree=cKDTree(midpoints+origin[np.newaxis,:])
        mp_tree =cKDTree(coords)
        # we need to span space with more boxes for small box size due to finite
        # representation of the contour points
        nbox=np.ceil(mx_dist/box_side)
        neighbour_boxes=mp_tree.query_ball_tree(
            box_tree,
            r=nbox*box_side
            )

        # neighbour_boxes=mp_tree.query_ball_point(midpoints,
        #                                          r=0.5*np.sqrt(2)*box_side
        #                                          )
        use_indices={ ii for nlist in neighbour_boxes for ii in nlist if nlist }
        boxes = [ makeSquare(midpoints[ii]+origin,box_side)
                    for ii in use_indices
                    ]
        # boxes = [ makeSquare(midpoint+origin,box_side) for midpoint in midpoints ]
        boxes = [ box for box in boxes if box.intersects(contourgon.exterior) ]
        side_lengths.append(box_side)
        num_boxes.append(len(boxes))

        # end=time.perf_counter()
        # print(f"{nn} divisions ({box_side=}) took {end-start} seconds")
        if nn==16:
            example_boxes=copy(boxes)

    # The gradient of the linear fit is the fractal dimension
    linear_fit = lambda x,m,c: m*x+c
    ln_num = np.log(num_boxes)
    ln_sides = -np.log(side_lengths)
    popt,pcov = scopt.curve_fit(linear_fit,
                                ln_sides,
                                ln_num
                                )
    # plt.plot(ln_sides,ln_num,
    #          marker='o',mfc='w',mec='k'
    #          )
    # plt.plot(ln_sides,linear_fit(ln_sides,*popt),'k--')
    # plt.show()
    # print(popt)
    return popt[0],example_boxes

def getContourCurvature(cells,coords,fname):
    """
        Find the local length scale
    """

    """
        Try to find the largest circle which fits in the colony
    """
    # find the centre of the colony and the rough radius
    ellipse=getColonyEllipse(coords)

    skel_name=fname.replace('rs_','skel_')
    print(f"{skel_name=}")
    try:
        skel_rs=np.load(skel_name)
    except Exception as e:
        print(e)

        # Select cells inside the centre of the colony
        dr=0.5
        density_field, density_grid = RodShapedBacterium.computeDensity(
            cells,
            dr=dr,
            fname=fname.replace('rs_','rho_')
            )

        # Threshold the density field to binarise the cell locationss
        eps=0.5
        density_field[density_field> eps]=1
        density_field[density_field<=eps]=0

        # Use medial axis to find the rough width of the colony at each location
        skel, distance = medial_axis(density_field,return_distance=True)
        dist_on_skel = distance*skel*dr
        if dist_on_skel.shape!=density_grid[0].shape:
            print(f"{dist_on_skel.shape=} {density_grid[0].shape=}")
            quit()
        # Filter anything too small and return the coordinates
        locs=np.where(dist_on_skel>1)
        dists=dist_on_skel[locs]
        X,Y=density_grid
        xs=X[locs]
        ys=Y[locs]
        skel_rs=np.array([xs,ys,dists]).T
        np.save(skel_name,skel_rs)

    dists=[]
    for rr in skel_rs:
        if Point(rr[:2]).within(ellipse):
            dists.append(rr[2])

    dists=np.array(dists)
    L2=dists.mean()
    print(f"{L2=}")

    return L2

def getCurvatureGridMethod(cells,coords,fname):
    try:
        rs=np.load(fname)
    except Exception as e:
        out_contour=computeColonyContourHoles(cells)

        smoothed_interiors=[]
        for interior in out_contour.interiors:
            p=Polygon(interior)
            if p.area>10:
                smoothed_interiors.append(interior)
        new_polygon = Polygon(out_contour.exterior.coords,holes=smoothed_interiors)

        left,down,right,up=out_contour.bounds
        origin=np.array([left,down])
        bounding_side_length=np.max([right-left,up-down])
        boundarygons=[ new_polygon, *[Polygon(n) for n in new_polygon.interiors] ]
        # boundarygons=[out_contour]

        nn=30
        xs,box_side = np.linspace(0,bounding_side_length,num=nn,
                                  endpoint=True,retstep=True
                                  )
        xx,yy = np.meshgrid(xs,xs)
        cx=xx.flatten()
        cy=yy.flatten()
        midpoints=np.array([mp for mp in zip(cx,cy)])
        pnts=[Point(midpoint+origin) for midpoint in midpoints]

        def distCost(p):
            pnt=Point(p)
            dst=lambda p: pnt.distance(p.exterior)
            closest_boundary=min(
                boundarygons,
                key=dst
                )
            return dst(closest_boundary)

        def invDistCost(p):
            return 1/distCost(p)

        ellipse=getColonyEllipse(coords)
        pnts=[ pp for pp in pnts if pp.within(ellipse) ]
        dsts=[ [pp.x,pp.y,distCost(pp)] for pp in pnts ]
        dsts=np.array(dsts)
        dsts=dsts[dsts[:,2]>1.0]
        print(f"{dsts.shape=}")
        ii=0
        count=0
        rs=[]
        fails=0
        for dd in dsts:
            pnt=dd[:2]
            res=scopt.minimize(invDistCost,x0=pnt)
            p=res.x
            if Point(p).within(ellipse):
                r=1/invDistCost(p)
                rs.append([*p,r])
                ii+=1
            else:
                fails+=1
        print(f"{count=} {fails=}")
        np.save(fname,rs)

    rs=np.array(rs)

    L1=rs[:,2].mean()
    print(f"{L1=}")
    return L1

def getContourCurvatureLegacy(cells,coords,fname):
    """
        Find the local length scale
    """

    """
        Try to find the largest circle which fits in the colony
    """
    skel_name=fname.replace('rs_','skel_')
    print(f"{skel_name=}")
    try:
        skel_rs=np.load(skel_name)
    except Exception as e:
        print(e)
        dr=0.5
        density_field, density_grid = RodShapedBacterium.computeDensity(
            cells,
            dr=dr
            )

        eps=0.5
        density_field[density_field> eps]=1
        density_field[density_field<=eps]=0
        skel, distance = medial_axis(density_field,return_distance=True)
        dist_on_skel = distance*skel*dr

        res=(distance==maximum_filter(distance,3,mode='constant',cval=0.0))
        locs=np.where((res==True)&(dist_on_skel>1))
        dists=dist_on_skel[locs]
        if dists.size==0:
            print("Filter too large, reducing to 2")
            res=(distance==maximum_filter(distance,2,mode='constant',cval=0.0))
            locs=np.where((res==True)&(dist_on_skel>1))
        dists=dist_on_skel[locs]
        X,Y=density_grid
        xs=X[locs]
        ys=Y[locs]
        skel_rs=np.array([xs,ys,dists]).T
        np.save(skel_name,skel_rs)

    # find the centre of the colony and the rough radius
    params=fitEllipse(coords)
    xc, yc, a, b, theta = params
    circ=Point([xc,yc]).buffer(1)
    ellipse=rotate(scale(circ,0.5*a,0.5*b),theta,use_radians=True)

    dists=[]
    for rr in skel_rs:
        if Point(rr[:2]).within(ellipse):
            dists.append(rr[2])
    # dists=skel_rs[:,2]
    # x,bins=np.histogram(dists,bins=20,density=True)
    # dr=np.diff(bins)[0]*0.5
    # L2=bins[np.argmax(x)]+dr
    # dists,_=filterOutliers(dists)
    dists=np.array(dists)
    L2=dists.mean()
    print(f"{L2=}")

    dr=0.5
    density_field, density_grid = RodShapedBacterium.computeDensity(
        [ cell for cell in cells if cell.getPolygon().intersects(ellipse) ],
        dr=dr
        )

    eps=0.5
    density_field[density_field> eps]=1
    density_field[density_field<=eps]=0
    skel, distance = medial_axis(density_field,return_distance=True)
    dist_on_skel = distance*skel*dr
    print(f"inside ellipse {np.mean(dist_on_skel[dist_on_skel>1])=}")
    quit()

    if False:
        try:
            rs=np.load(fname)
        except Exception as e:
            out_contour=computeColonyContourHoles(cells)

            smoothed_interiors=[]
            for interior in out_contour.interiors:
                p=Polygon(interior)
                if p.area>25:
                    smoothed_interiors.append(interior)
            new_polygon = Polygon(out_contour.exterior.coords,holes=smoothed_interiors)

            left,down,right,up=out_contour.bounds
            origin=np.array([left,down])
            bounding_side_length=np.max([right-left,up-down])
            boundarygons=[ new_polygon, *[Polygon(n) for n in new_polygon.interiors] ]
            # boundarygons=[out_contour]

            def invDistCost(p):
                pnt=Point(p)
                key=lambda p: pnt.distance(p.exterior)
                closest_boundary=min(
                    boundarygons,
                    key=key
                    )
                return 1/key(closest_boundary)

            ii=0
            count=0
            rs=[]
            fails=0
            while ii<30:
                count+=1
                pnt=origin+np.array([right-left,up-down])*np.random.random(2)
                if Point(pnt).within(out_contour):
                    res=scopt.minimize(invDistCost,x0=pnt)
                    p=res.x
                    if Point(p).within(out_contour):
                        r=1/invDistCost(p)
                        rs.append([*p,r])
                        ii+=1
                    else:
                        fails+=1
            print(f"{count=} {fails=}")
            np.save(fname,rs)

        rs=np.array(rs)
        # x,bins=np.histogram(rs[:,2],bins=20,density=True)
        # dr=np.diff(bins)[0]*0.5
        # L1=bins[np.argmax(x)]+dr
        # rds=rs[rs[:,2]>6]
        # rds,_=filterOutliers(rds[:,2])
        L1=rs[:,2].max()
        print(f"{L1=}")
        return L1,L2,rs,skel_rs
    return L2

# def findCentreLine(cells):
#     out_contour=computeColonyContourHoles(cells).exterior
#     bounds=out_contour.bounds
#
#     quit()

def getColonyWavelength(cells_,nbins=10):
    """
    Try to extract an indication of the global wavelength

    Plan:
        Bin slices in some x_range by their COM
        Fit a Bezier curve to this
        Extract something like wavelength from the Bezier
    """

    print("Place a point at set x with maximum distance from edges")
    print("This will then be fit with something")
    quit()
    # Sort the cells by their x coordinate
    cells_.sort(key=lambda c: c.rcm[0])
    x_min=cells_[0].rcm[0]
    x_max=cells_[-1].rcm[0]
    edges=np.linspace(x_min,x_max+0.01,nbins+1)
    inds=np.digitize([c.rcm[0] for c in cells_],bins=edges)
    # for ii in range(len(cells_)):
    #     print(edges[inds[ii]-1], "<=", cells_[ii].rcm, "<", edges[inds[ii]])

    rcms=[]
    thetas=[]
    for ee in range(1,len(edges)):
        cell_inds=np.where(inds==ee)[0]
        rcms.append(np.array([cells_[ii].rcm for ii in cell_inds]).mean(axis=0) )
        thetas.append(np.array([cells_[ii].theta for ii in cell_inds]).mean() )
    rcms=np.array(rcms)
    thetas=np.array(thetas)
    thetas-=thetas.mean()
    thetas=np.convolve(thetas,np.full(10,1/10))

    plt.plot(thetas)
    plt.show()

    pi=scipy.signal.find_peaks(np.abs(thetas),height=0.5*np.max(thetas),width=5,
                               distance=len(thetas)//5)
    print(pi)
    pi=pi[0]
    rough_lam=2*np.diff(pi).mean()
    print(f"{rough_lam=}")
    xs=rcms[:,0]
    # ys=rcms[:,1]
    # f=lambda x,a,x0,c,lam: a*np.sin(2*np.pi*(x-x0)/lam)+c
    # popt,pcov=scipy.optimize.curve_fit(f,xs,thetas,sigma=np.full(len(thetas),0.1))
    # print(f"{popt=}")
    # fig,ax=plt.subplots(1,1)
    # addAllCellsToPlot(cells_,ax,ax_rng=x_max-x_min,
    #                   alpha=0.8,show_id=False)
    # ax.scatter(xs,ys,s=20,fc='w',ec='k',alpha=0.7)
    # ax.plot(xs,f(xs,*popt),'b--',lw=2)
    # plt.show()

    ftheta=scipy.fft.rfft(thetas)
    # ixmax=np.argpartition(np.abs(ftheta),-1)[-1:]
    ixmax=np.argmax(np.abs(ftheta))
    print(f"fourier lam: {len(thetas)/ixmax}")

    f2=lambda xs,a: np.sum(
        np.array([a])[:,None]*np.array(
            [np.sin(xs*2*np.pi*ii/len(thetas)) for ii in [ixmax]]
            ),
        axis=0
        )
    popt,pcov=scipy.optimize.curve_fit(f2,xs,thetas)
    print(f"{popt=}")
    fig,ax=plt.subplots(1,2)
    ax[0].plot(xs,thetas)
    ax[0].scatter(xs[pi],thetas[pi])
    # ax[0].plot(xs,f2(xs,*popt))
    ax[0].plot(xs,np.sin(2*np.pi*xs/rough_lam),':')

    ax[0].plot(xs,np.sin(xs*2*np.pi*ixmax/len(thetas)),'--')
    ax[1].plot(2*np.pi*np.arange(1+len(thetas)//2)/len(thetas),np.abs(ftheta))
    # ax.plot(xs,f(xs,*popt),':')
    plt.show()
    quit()
    return popt,thetas,rcms

def getCorrelationLengths(cells,dr,rho_file,Q_file,C_file):
    """
    cells: cell list
    dr: fineness of the grid
    data_file: expected to be a binary format to try to load data if prev. calc'd

    q is the charge field, Q is the tensor and Q_grid is the domain in XX,YY
    """

    dir,fname=os.path.split(C_file)
    fname=fname.replace('Correlation','distance_list',1)
    d_grid_file=f"{dir}/{fname}"
    try:
        C_grid=np.load(C_file)
        d_grid=np.load(d_grid_file)
    except Exception as e:
        rho, grid = RodShapedBacterium.computeDensity(cells,dr=dr,fname=rho_file)
        q_field, Q_field, Q_grid = RodShapedBacterium.computeChargeDensity(
            cells,
            dr=dr,
            fname=Q_file
            )
        X,Y=Q_grid
        # d_grid=np.meshgrid(np.arange(-(nx//2),nx//2+1),np.arange(-(ny//2),ny//2+1))
        nn=1+int(50/dr)
        xv,yv=np.meshgrid(np.arange(-(nn),nn+1),np.arange(-(nn),nn+1))
        d_grid=np.array([xv.flatten(),yv.flatten()],dtype=int).T
        C_grid=np.zeros((d_grid.shape[0],2),dtype=float)

        set_num_threads(12)
        findCorrelationLengths(rho,Q_field,dr,d_grid,C_grid)
        print(findCorrelationLengths.parallel_diagnostics(level=4))
        C_grid=C_grid/np.sum(rho)
        np.save(C_file,C_grid)
        np.save(d_grid_file,d_grid)

    ds=dr*np.linalg.norm(d_grid,axis=1)
    c_par=C_grid[:,0]
    c_perp=C_grid[:,1]

    cell_diam=1
    print(f"{ds.max()=}")
    num_bins=int(np.ceil(ds.max()/cell_diam))
    edges=np.arange(0,(num_bins+1)*cell_diam,cell_diam)
    avgs=np.zeros((num_bins,2))
    print(edges)
    for ii in range(1,num_bins+1):
        idx=((ii-1)*cell_diam<ds) & (ds<ii*cell_diam)
        avgs[ii-1,0]=c_par[idx].mean()
        avgs[ii-1,1]=c_perp[idx].mean()
    print(avgs)

    plt.plot(edges[:-1],avgs[:,0],label='par')
    plt.plot(edges[:-1],avgs[:,1]/avgs[0,0],label='perp')
    plt.show()

    l_par =scipy.integrate.trapezoid(avgs[:,0],dx=cell_diam)/avgs[0,1]
    l_perp=scipy.integrate.trapezoid(avgs[:,1],dx=cell_diam)/avgs[0,0]

    print(f"{l_par=} {l_perp=}")
    # print(scipy.integrate.simpson(c_par)/c_perp[0])
    # plt.close('all')
    # plt.imshow(rho)
    # plt.show()
    quit()
    return l_par,l_perp

def getPairCorrelationData(cells):
    N=len(cells)
    rs=np.zeros((N*(N-1)//2),dtype=float)-1 # modulus distance
    phis=np.zeros((N*(N-1)//2),dtype=float) # order param
    cell_dat = RodShapedBacterium.toNative(cells)
    findPairCorrelationNumba(cell_dat,rs,phis)
    return rs,phis

@njit
def findPairCorrelationNumba(cell_dat,rs,phis):
    # for cc,pnt in enumerate(result):
    # norm_factor=0.0
    N=cell_dat.shape[0]
    count=0
    for ii in range(N):
        rcm1=cell_dat[ii,0:2]
        tang1=cell_dat[ii,2:4]
        length1=cell_dat[ii,4]

        for jj in range(ii):
            rcm2=cell_dat[jj,0:2]
            tang2=cell_dat[jj,2:4]
            length2=cell_dat[jj,4]

            rs[count]=np.linalg.norm(rcm1-rcm2)
            phis[count]=2*np.dot(tang1,tang2)**2-1
            # *length1*length2
            # norm_factor+=length1*length2
            count+=1
    # phis/=norm_factor

def getWeightedPairCorrelationData(cells,width=1):

    N=len(cells)
    ax_lim = np.max(VisInfectedBiofilm.getAxisLimits(cells)) + 5
    # phis=np.zeros((N*(N-1)//2),dtype=float) # order param
    # weights=np.zeros((N*(N-1)//2),dtype=float) # average weights
    # bin_edges=np.arange(0,ax_lim+width,width)

    # RodShapedBacterium.getNeighbourLists(cells,r=ax_lim,find_contacts=False)
    # rs=np.arange(0,R+diameter,20*diameter)
    # rs=np.array([0,*np.logspace(0,np.log2(ax_lim-4),base=2,num=8)])
    rs=np.array([*np.arange(0,5+width,width),
        *np.logspace(np.log2(5+width),np.log2(ax_lim-4),base=2,num=10)
        ])
    num_bins=len(rs)-1
    print(f"{num_bins=}")
    print(rs)


    colony_contour=computeColonyContour(cells)
    ellipse=getColonyEllipse(colony_contour)

    ref_cells=[ cc for cc in cells if ellipse.contains(cc.getPolygon()) ]
    Nref=len(ref_cells)
    c_ls=np.zeros((Nref,num_bins))
    avgs=np.zeros(num_bins)

    ngons=[ [nn,nn.getPolygon()] for nn in cells ]
    refgons=[ [nn,nn.getPolygon()] for nn in ref_cells ]

    for ii in range(Nref):
        cc_data=refgons[ii]
        cc=cc_data[0]
        if ellipse.contains(cc_data[1])==False:
            continue
        centre=cc.rcm[:-1]

        ds=np.array([ np.linalg.norm(x[0].rcm[:-1]-centre) for x in ngons ])
        indxs=np.argsort(ds)
        ds=ds[indxs]
        ngons=[ ngons[jj] for jj in indxs ]
        for bb in range(1,num_bins+1):
            # make an annulus of width diam and outer radius r
            # bin cells into this annulus
            # find the average dot product weighted by area in the annulus
            r2=rs[bb]
            r1=r2-width
            ring=makeAnnulus(centre,r1=r1,r2=r2)

            lower_lim=np.argmax(ds>=(r1-3))
            upper_lim=np.argmax(ds>=(r2+3))
            # print(f"{lower_lim=} {upper_lim=} {r1=} {r2=}")
            binned_ngons=[ ngons[jj] for jj in range(lower_lim,upper_lim+1)
                           if ngons[jj][1].intersects(ring) ]
            if not binned_ngons:
                continue

            # fig,ax=plt.subplots()
            # # RodShapedBacterium.sus_vis_radius_factor=0.7
            # fcp.addAllCellsToPlot(cells,ax,ax_lim/2,alpha=0.3)
            # ax.add_patch(PolygonPatch(ring,fc="None"))
            # for gg in binned_ngons:
            #     # ax.plot(*gg[1].exterior.coords.xy)
            #     try:
            #         mp=gg[1].intersection(ring)
            #         ax.add_patch(PolygonPatch(mp,fc=RodShapedBacterium.colour_fn(gg[0].theta)))
            #     except Exception as e:
            #         print(e)
            #         mp_intersection=gg[1].intersection(ring)
            #         print(f"{mp_intersection.area=} {[a.area for a in mp_intersection.geoms]}")
            #         for mp in mp_intersection.geoms:
            #             ax.add_patch(PolygonPatch(mp,fc=RodShapedBacterium.colour_fn(gg[0].theta)))
            # ax.add_patch(PolygonPatch(cc.getPolygon(),fc='None',hatch='//',ec='k',ls='--'))
            #
            # ax.axis("scaled")
            # ax.axis([*np.array(ring.bounds)[[0,2,1,3]]])
            # fig.savefig(f"/home/rory/example_correlation/{ii}_rings_{bb}_2.png",
            #             format="png",dpi=200,bbox_inches='tight',
            #             transparent=True)
            # plt.show()
            # plt.close()

            intersections=[bng[1].intersection(ring).area
                            for bng in binned_ngons]
            total_area=sum(intersections)

            def f(a,b,intersect):
                dp=np.dot(a.ori,b[0].ori)
                return (2*dp**2-1)*intersect

            # make this a nematic dot product
            c_ls[ii,bb-1]=sum([f(cc,binned_ngons[kk],intersections[kk])
                              for kk in range(len(binned_ngons))])/total_area
        # quit()

    c_ls=c_ls[c_ls[:,0]>0].mean(axis=0)
    # c_ls/=np.max(c_ls)
    rs=rs[1:]-0.5*width
    # plt.plot(rs,c_ls)

    # scal=lambda x,a,b,xi: a*np.exp(-x/xi)+b
    # popt,pcov=scipy.optimize.curve_fit(scal,rs,c_ls,p0=(1,0,30))
    # plt.plot(rs,scal(rs,*popt),'--')
    # plt.show()
    # quit()

    # N=len(cells)
    # ax_lim = np.max(VisInfectedBiofilm.getAxisLimits(cells)) + 5
    # phis=np.zeros((N*(N-1)//2),dtype=float) # order param
    # weights=np.zeros((N*(N-1)//2),dtype=float) # average weights
    # bin_edges=np.arange(0,ax_lim+width,width)
    #
    # num_bins=len(bin_edges)-1
    # print(f"{num_bins=}")
    # c_ls=np.zeros((N,num_bins))
    # avgs=np.zeros(num_bins)
    # weights=np.zeros(num_bins)
    #
    # cell_gons=[ nn.getPolygon() for nn in cells ]
    # rings=[ makeAnnulus([0.0,0.0],r1=bin_edges[bb-1],r2=bin_edges[bb])
    #         for bb in range(1,num_bins) ]
    # # affinity.translate
    # # quit()
    # trials=10
    # for _ in tqdm(range(trials)):
    #     ii = np.random.randint(low=0,high=N-1,dtype=int)
    #     cc=cells[ii]
    #     for jj in range(N):
    #         if ii==jj:
    #             continue
    #         nn=cells[jj]
    #         centre=cc.rcm[:-1]
    #         ngon=cell_gons[jj]
    #         r_lims=np.array([ nn.rcm[:-1]+nn.ori[:-1]*nn.length*0.5,
    #                           nn.rcm[:-1]-nn.ori[:-1]*nn.length*0.5 ])
    #         r_lims=np.sort(np.linalg.norm(r_lims-centre,axis=1))
    #
    #         ind=np.digitize(r_lims,bins=bin_edges)
    #
    #         for bb in range(ind[0],ind[1]+1):
    #             # make an annulus of width diam and outer radius r
    #             # bin cells into this annulus
    #             # find the average dot product weighted by area in the annulus
    #             # r2=bin_edges[bb]
    #             # r1=bin_edges[bb-1]
    #             # print(f"{r1=} {r2=}")
    #             # ring=makeAnnulus(centre,r1=r1,r2=r2)
    #             ring=translate(rings[bb-1],xoff=centre[0],yoff=centre[1])
    #
    #             area_poly=ring.intersection(ngon)
    #
    #             area=area_poly.area
    #             weights[bb-1]+=area
    #             avgs[bb-1]+=area*(2*np.dot(cc.ori[:-1],nn.ori[:-1])**2-1.0)
    #
    #             # fig,ax=plt.subplots()
    #             # wedge_patch = PolygonPatch(ring,fc="None")
    #             # ax.add_patch(wedge_patch)
    #             # try:
    #             #     ax.add_patch(PolygonPatch(area_poly))
    #             # except Exception as e:
    #             #     for mp in area_poly.geoms:
    #             #         ax.add_patch(PolygonPatch(mp))
    #             # ax.axis('scaled')
    #             # ax.axis([-100,100,-100,100])
    #             # plt.show()
    #             # quit()
    #             #
    #             # # make this a nematic dot product
    #             # c_ls[ii,bb-1]=sum([f(cc,bng) for bng in binned_ngons])/total_area
    # c_ls=avgs/weights
    # rs=bin_edges[:-1]
    #
    # inds=np.argwhere(~np.isnan(c_ls))
    # c_ls=c_ls[inds]
    # rs=rs[inds]
    #
    # print(f"{c_ls=}")
    #
    # scal=lambda x,xi: np.exp(-x/xi)
    # popt,pcov=scipy.optimize.curve_fit(scal,rs,c_ls,p0=(1))
    # plt.plot(rs,scal(rs,*popt),'--')
    # plt.show()
    # quit()
    #
    #
    # return rs,phis,weights
    return rs,c_ls

@njit
def findPairCorrelationNumba(cell_dat,rs,phis):
    # for cc,pnt in enumerate(result):
    # norm_factor=0.0
    N=cell_dat.shape[0]
    count=0
    for ii in range(N):
        rcm1=cell_dat[ii,0:2]
        tang1=cell_dat[ii,2:4]
        length1=cell_dat[ii,4]

        for jj in range(ii):
            rcm2=cell_dat[jj,0:2]
            tang2=cell_dat[jj,2:4]
            length2=cell_dat[jj,4]

            rs[count]=np.linalg.norm(rcm1-rcm2)
            phis[count]=np.dot(tang1,tang2)**2-0.5
            # *length1*length2
            # norm_factor+=length1*length2
            count+=1
    # phis/=norm_factor


@njit
def calcCorrIntegrandNumba(rho,Q,Qp,dhat):
    Qxx=Q[0,0]
    Qxy=Q[0,1]
    Qpxx=Qp[0,0]
    Qpxy=Qp[0,1]
    g=Qxx*Qpxx+Qxy*Qpxy
    parallel=1+dhat[0]*dhat[0]*Qxx+dhat[0]*dhat[1]*Qxy
    perp=1+dhat[0]*dhat[0]*Qxx-dhat[0]*dhat[1]*Qxy
    return np.array([rho*g*parallel,rho*g*perp])

@njit(parallel=True)
def findCorrelationLengths(rho,Q_field,dr,d_grid,C_grid):
    # si,sj=d_grid.shape
    nx,ny=rho.shape
    out_freq=d_grid.shape[0]//20
    for cc in range(d_grid.shape[0]):
        if cc%out_freq==0:
            print("progress: ",100*cc/d_grid.shape[0])
        d=d_grid[cc]

        # Get the shift
        di,dj=d
        if di==0 and dj==0:
            dhat=np.array([0.0,0.0])
        else:
            dhat=d/np.sqrt(np.sum(d**2))

        # Integrate over the grid
        total=np.array([0.0,0.0])
        for ii in prange(nx):
            for jj in range(ny):
                xprime=ii+di
                yprime=jj+dj
                if xprime<nx and yprime<ny and xprime>=0 and yprime>=0:
                    total+=calcCorrIntegrandNumba(
                        rho[ii,jj],
                        Q_field[ii,jj],
                        Q_field[xprime,yprime],
                        dhat
                        )
        C_grid[cc]=total

def getAreaFraction(cells):

    # findCentreLine(cells)
    #
    # quit()
    out_contour=computeColonyContourHoles(cells)
    # smoothed_interiors=[]
    # for interior in out_contour.interiors:
    #     p=Polygon(interior)
    #     if p.area>5+np.pi:
    #         smoothed_interiors.append(interior)
    # hole_area=np.sum([Polygon(p).area for p in smoothed_interiors])
    # new_polygon = Polygon(out_contour.exterior.coords,holes=smoothed_interiors)

    area_frac=out_contour.area/Polygon(out_contour.exterior).area
    return area_frac

def checkBuckled(cells):

    area_frac=getAreaFraction()

    return area_frac<0.99

def getCriticalLength():
    print("getCriticalLength quittiing")
    quit()
    #===
    """
    Set up skeletonize
    """
    dr=0.5
    density_field, density_grid = RodShapedBacterium.computeDensity(
        cells,
        dr=dr
        )
    density_field[density_field>0.05]=1
    density_field[density_field<0.05]=0
    skeleton=skeletonize(density_field)
    fig,ax=plt.subplots()
    cs=ax.contourf(skeleton,
                   cmap=cm.PuBu_r)
    cbar=fig.colorbar(cs)
    plt.show()
    quit()

    #===
    """
    Try to find the local wavelength distribution
    """
    dr=0.5
    density_field, density_grid = RodShapedBacterium.computeDensity(
        cells,
        dr=dr
        )
    N,M=density_field.shape

    fig,ax=plt.subplots()
    # for ii in range(0,N,N//20):
    #     q=np.fft.fft(density_field[ii])
    #     ax.plot()
    X,Y=density_grid
    qs=np.abs(np.fft.rfft2(density_field))
    N,M=qs.shape
    X=np.arange(N)*2*np.pi/(N*dr)
    Y=np.arange(M)*2*np.pi/(M*dr)
    cs=ax.contourf(X,Y,qs.T,locator=ticker.LogLocator(),cmap=cm.PuBu_r)
    cbar=fig.colorbar(cs)
    ax.set_ylabel("$q_y$ \si{\per\micro\meter}")
    ax.set_xlabel("$q_x$ \si{\per\micro\meter}")
    plt.show()
    quit()




    #===
    # Aim to find the distribution of circles we can fit into the contour
    contourgon = Polygon([(v[0],v[1]) for v in coords])
    pnt = np.random()
    quit()
    coords = computeColonyContour(cells)
    xs=coords[:,0]
    ys=coords[:,1]


    # theta=np.linspace(0,2*np.pi,num=10000)
    # xs=np.cos(theta)
    # ys=np.sin(theta)
    plt.plot(xs,ys,label="contour",alpha=0.8)

    # n=100
    # v=np.full(n,1/n)
    # smooth = lambda x : np.convolve(x,v,mode='valid')
    # xs = smooth(xs)
    # ys = smooth(ys)
    # plt.plot(xs,ys,'o')

    # Interpolate
    tck,u=splprep([xs,ys],per=True,k=3,s=1000)
    # unew=np.linspace(u[0],u[-1],len(u))
    unew=u

    new_points=splev(unew,tck)
    plt.plot(new_points[0],new_points[1],label=f"splprep",alpha=0.8)

    # Evaluate curvature
    xp,yp=splev(unew,tck,der=1)
    xpp,ypp=splev(unew,tck,der=2)
    K=np.abs(xp*ypp-yp*xpp)/np.power(xp**2+yp**2,3/2)

    ins=np.logical_and(0.5e-2<K,K<0.1)
    plt.scatter(new_points[0][ins],new_points[1][ins],c=1/K[ins],
                cmap=cm.RdYlGn,
                s=100,edgecolor="None",label="$K$",alpha=0.8)
    plt.legend()
    plt.show()
    plt.plot(unew,K)
    plt.show()
    quit()
