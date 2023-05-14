# Standard modules
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import pandas as pd
import os

# Third party
from shapely.geometry import Polygon

# Custom modules
from RodShapedBacteria import RodShapedBacterium
from findCriticalLengthProb import getCells
import utilities as ut
import fastCellPlotting as fcp
from DistributionFunctions import computeColonyContour,getColonyDensity
from generalPlotting import addNematicDirector

save_dir="../GeneratedOutput/AnalysisResults/ScalingOrientation"
# bdir="../GeneratedOutput/SimOutput/testStiffness/"
# bdir+="radius_30_width_9_linking_prob_0_kappa_0_B_0/"
# data_dir="../GeneratedOutput/SimOutput/OrientationalScalingExpCheckNoRelax/"
data_dir="../GeneratedOutput/SimOutput/testStiffness/"

add_boundary=False

for dd in glob(f"{data_dir}/*"):
    for channel in ['annulus','straight']:
        RodShapedBacterium.sus_vis_radius_factor=1
        file_initial=f"{dd}/biofilm_{channel}_channel_initial_00000.dat"
        file_relaxed=f"{dd}/biofilm_{channel}_channel_relaxed_00000.dat"
        hps=pd.read_csv(f"{dd}/LcDist.log",delimiter='\t')
        fname=dd.split('/')[-1]

        for ii,file in enumerate([file_initial,file_relaxed]):
            cells=getCells(file)
            fig,ax=plt.subplots(1,1)

            for cell in cells:
                cell.colour=RodShapedBacterium.colour_fn(cell.theta)
                # if cell.cell_id==38 or cell.cell_id==34:
                #     cell.addElementToPlot(ax,colour=None,projection="xy",ax_rng=80,show_id=True)
            gon=computeColonyContour(cells,add_links=True,ret_gon=True)
            print(f"{gon.area=}")

            if add_boundary:
                ax.plot(*gon.exterior.xy,color='k',lw=3,alpha=0.4)
            for geom in gon.interiors:
                if geom.length>=0.1*gon.exterior.length:
                    inner=geom
                    if add_boundary:
                        ax.plot(*geom.xy,color='k',lw=3,alpha=0.4)
            try:
                outer=Polygon(np.array([(z[0],z[1]) for z in zip(*gon.exterior.xy)]))
                inner=Polygon(np.array([(z[0],z[1]) for z in zip(*inner.xy)]))
                annulus=outer.difference(inner)
            except Exception as e:
                print(e)
                annulus=Polygon(np.array([(z[0],z[1]) for z in zip(*gon.exterior.xy)]))

            print(f"{annulus.area=}")
            density=sum([cell.getArea() for cell in cells])/annulus.area

            minx, miny, maxx, maxy=gon.bounds
            fcp.addAllCellsToPlot(cells,ax,ax_rng=maxx-minx,show_id=True)
            fcp.addChainLinks(cells,ax,ax_rng=maxx-minx,show_id=True)
            ax.axis('scaled')
            ax.set_xlim([minx,maxx])
            ax.set_ylim([miny,maxy])

            mod=''
            if ii==0:
                mod='initial'
            elif ii==1:
                mod='relaxed'
            else:
                quit()
            ax.text(0,0,fr"$\phi={density:3.3f}$ {mod}",ha='center',va='center')
            ax.axis('off')
            q_name=f"{save_dir}/{fname}_Q_{channel}_{mod}_vis.npy"
            addNematicDirector(ax,cells,q_name)
            fig.savefig(f"{save_dir}/{fname}_{channel}_{mod}_vis.png",
                        dpi=300,
                        transparent=True,
                        bbox_inches='tight'
                        )
            fig.savefig(f"{save_dir}/{fname}_{channel}_{mod}_vis.pdf",
                        transparent=True,
                        bbox_inches='tight'
                        )
        plt.show()
