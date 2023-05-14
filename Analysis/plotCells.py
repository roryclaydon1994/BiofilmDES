# Standard modules
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import pandas as pd
import os
import sys

# Third party
from shapely.geometry import Polygon

# Custom modules
from RodShapedBacteria import RodShapedBacterium
from findCriticalLengthProb import getCells
import utilities as ut
import fastCellPlotting as fcp
from DistributionFunctions import computeColonyContour,getColonyDensity
from generalPlotting import addNematicDirector

ut.setMPL()
save_dir="/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/figures/directors/"

def plotCells(file,add_director=True):
    RodShapedBacterium.sus_vis_radius_factor=0.7
    streamplot=(True & add_director)

    fig,ax=plt.subplots(1,1,figsize=ut.getFigsize(cols=2),
                        constrained_layout=True,
                        facecolor='k')
    cells=getCells(file)
    for cell in cells:
        cell.colour=RodShapedBacterium.colour_fn(cell.theta)

    dat=pd.read_csv(file,sep='\t')
    maxx=dat['pos_x'].max()+3
    minx=dat['pos_x'].min()-3
    maxy=dat['pos_y'].max()+3
    miny=dat['pos_y'].min()-3
    fcp.addAllCellsToPlot(cells,ax,ax_rng=maxx-minx,show_id=True,ec='w')
    fcp.addChainLinks(cells,ax,ax_rng=maxx-minx,show_id=True,colour='w')
    if add_director:
        q_name=f"{save_dir}/"
        q_name+=f"{file.replace('/','_')}_Q.npy"
        addNematicDirector(ax,cells,q_name,streamplot=streamplot,dr=3)
    ax.axis('scaled')
    ax.axis('off')
    ax.set_xlim([minx,maxx])
    ax.set_ylim([miny,maxy])
    stem=f"{save_dir}/plotCells_{file.replace('/','_')}"
    fig.savefig(f"{stem}_{streamplot=}.png",
                dpi=300,
                transparent=False,
                bbox_inches='tight'
                )
    fig.savefig(f"{stem}_{streamplot=}.pdf",
                transparent=False,
                bbox_inches='tight'
                )
    plt.show()
if __name__ == "__main__":
    if int(sys.argv[2])==0:
        add_director=False
    else:
        add_director=True
    plotCells(sys.argv[1],add_director=add_director)
