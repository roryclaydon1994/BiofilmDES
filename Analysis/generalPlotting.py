##

# Standard modules
import matplotlib.pyplot as plt
import numpy as np
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Third party modules
from numba import jit,njit,prange,set_num_threads

# Custom modules
import fastCellPlotting as fcp
from RodShapedBacteria import RodShapedBacterium

def plotUnstructuredMesh(cells,field,coords,fname):

    RodShapedBacterium.sus_vis_radius_factor=0.7

    xs=coords[:,0]
    ys=coords[:,1]

    # xi = np.linspace(xs.min()-5,xs.max()+5, 200)
    # yi = np.linspace(ys.min()-5,ys.max()+5, 200)

    ax_lim=coords.max()-coords.min()

    # From https://matplotlib.org/stable/gallery/images_contours_and_fields/irregulardatagrid.html
    # Linearly interpolate the data (x, y) on a grid defined by (xi, yi).
    # ax.tricontour(x, y, field, levels=14, linewidths=0.5, colors='k')

    fig,ax=plt.subplots(1,1,figsize=[5,5])
    ax.axis('off')
    ax.axis('scaled')
    ax.axis([xs.min(),xs.max(),ys.min(),ys.max()])

    CS = ax.tricontourf(xs, ys, field, levels=14, cmap="PRGn")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar=fig.colorbar(CS,cax=cax)

    cbar.set_label("$S$",size=14)
    cbar.ax.tick_params(labelsize=12)
    for cell in cells:
        cell.colour="None"

    fcp.addAllCellsToPlot(
        cells=cells,
        ax=ax,
        ax_rng=ax_lim,
        show_id=False
        )

    fig.savefig(fname,
                format=os.path.splitext(fname)[1].replace('.', ''),
                bbox_inches='tight',
                transparent=True,
                dpi=300)

@njit(parallel=True)
def getNs(nematic_director,Q,ny,nx):
    for yi in prange(ny):
        for xi in prange(nx):
            w,v=np.linalg.eig(Q[yi,xi])
            pos_index=np.argsort(w)[-1]
            if w[pos_index]<1e-1:
                continue
            n=v[:,pos_index]
            nematic_director[yi,xi]=n/np.linalg.norm(n)

def addNematicDirector(ax,cells,q_name,dr=1,streamplot=False,start_point=None):
    q_name=q_name.replace('Q',f'Q_{dr=}')
    q,Q,grid=RodShapedBacterium.computeChargeDensity(
        cells,dr=dr,
        fname=q_name
        )
    ny,nx=grid[0].shape
    nematic_director=np.zeros( (ny,nx,2) )
    getNs(nematic_director,Q,ny,nx)
    U=nematic_director[:,:,0]
    V=nematic_director[:,:,1]

    if streamplot:
        strm=ax.streamplot(grid[0],grid[1],U,V,broken_streamlines=True,
                           integration_direction='both',
                           linewidth=1,
                           color='w',
                           arrowsize=0,
                           # minlength=4,
                           # maxlength=np.inf,
                           density=2.5)
        strm.lines.set_colors(np.array([1,1,1,1]))
        print(strm.lines.get_colors())
        return strm
    else:
        mask = np.logical_or(V!=0,U!=0)
        quvr=ax.quiver(grid[0][mask],grid[1][mask],U[mask],V[mask],pivot='middle',
                       angles='xy',scale_units='xy',scale=None,
                       headaxislength=0,headlength=0,headwidth=1
                       )
        return quvr
