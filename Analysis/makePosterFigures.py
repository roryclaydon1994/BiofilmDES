# Standard modules
import argparse
import numpy as np
import pandas as pd
import scipy
import scipy.stats
from glob import glob
from pprint import pprint
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import os
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
fontprops = fm.FontProperties(size=36)

# Custom modules
from visualise_biofilm import setFigSize
from ChainingRodShapedBacteria import ChainingRodShapedBacterium
from RodShapedBacteria import RodShapedBacterium, getNematicAngle
from FindDefects import intialiseElementsFromData
from visualise_infected_biofilm import VisInfectedBiofilm
from findCriticalLengthProb import getCells
import utilities as ut
import fastCellPlotting as fcp

ut.setMPL()

def makeColourWheel(colour_fn):
    # plt.close()
    fig,ax = plt.subplots(1,1,figsize=[10,10])
    angles = np.linspace(-np.pi, np.pi, 360)
    ax = plt.subplot(1, 1, 1, polar=True)
    for sa,ea in zip(angles[:-1],angles[1:]):
        ax.fill_between([sa,ea],0,1,color=colour_fn(sa))
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_axis_off()
    base_dir="/home/rory/PhD_sync/BiofilmPhageDES/poster"
    fig.savefig(f"{base_dir}/colour_wheel.pdf",
                format="pdf",
                bbox_inches='tight',
                dpi=dpi,
                transparent=True)
    # plt.show()

def addCellsPoster(file,ax,shift,ax_rng):
    RodShapedBacterium.sus_vis_radius_factor=0.7
    cells=getCells(file)
    for cell in cells:
        cell.colour=RodShapedBacterium.colour_fn(cell.theta)
        cell.rcm[0]+=shift[0]
        cell.rcm[1]+=shift[1]


    dat=pd.read_csv(file,sep='\t')

    fcp.addAllCellsToPlot(cells,ax,ax_rng=ax_rng,show_id=True)
    fcp.addChainLinks(cells,ax,ax_rng=ax_rng,show_id=True,colour='w')

def createImage(rn,files,bounds,ax_rng,format='png'):
    min_min_x,max_max_x,min_min_y,max_max_y=bounds
    fig=plt.figure(facecolor='k')
    ax=plt.axes()
    for idx,file in enumerate(files):
        # file=glob(f"{base_dir}/run{ii}/repeat3/final*")[0]
        r=idx%2
        c=idx//2
        shift=np.array([r*500,c*350])
        addCellsPoster(file,ax,shift,ax_rng)
        ax.set_xlim([min_min_x,max_max_x])
        ax.set_ylim([min_min_y,max_max_y])
        # file=glob(f"{base_dir}/run{ii}/repeat1/biofilm_00138.dat")[0]
        # print(f"Loading from {file}")
        # plot_dict[ii]=getCells(file)
        # clim=np.max(VisInfectedBiofilm.getAxisLimits(plot_dict[ii]))+5
    ax.axis('scaled')
    ax.axis('off')
    # ax.set_facecolor('k')
    ax.set_xlim([min_min_x,max_max_x])
    ax.set_ylim([min_min_y,max_max_y])
    base_dir=f"/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/NBIC/animation{run_nums}"
    # save_name=f"{base_dir}/NBIC_plot.pdf"
    # fig.savefig(save_name,
    #             format="pdf",
    #             bbox_inches='tight',
    #             dpi=dpi,
    #             # transparent=True
    #             )
    save_name=f"{base_dir}/biofilm_{rn:05d}.{format}"
    fig.savefig(save_name,
                format=format,
                bbox_inches='tight',
                dpi=400,
                # transparent=True
                )
    # save_name=f"{base_dir}/NBIC_plot.png"
    # fig.savefig(save_name,
    #             format="png",
    #             bbox_inches='tight',
    #             dpi=dpi,
    #             # transparent=True
    #             )
    # plt.show()
    plt.close()

def createSingleImage(run_num,outnum):
    base_dir="/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/"
    base_dir+="SimOutput/ChainingPhaseDiagramFinal"
    ax_lims=0

    file=f"{base_dir}/run{run_num}/repeat3/biofilm_{outnum:05d}.dat"
    dat=pd.read_csv(file,sep='\t')

    max_max_x=dat['pos_x'].max()+2.5
    min_min_x=dat['pos_x'].min()-2.5
    max_max_y=dat['pos_y'].max()+2.5
    min_min_y=dat['pos_y'].min()-2.5

    ax_rng=max_max_x-min_min_x
    print(f"{ax_rng=}",f"{max_max_x=}",f"{min_min_x=}",f"{max_max_y=}",f"{min_min_y=}")
    files=[file]
    bounds=(min_min_x,max_max_x,min_min_y,max_max_y)
    createImage(run_num,files,bounds,ax_rng,format='pdf')

def makeComparisonMorphology():
    base_dir="/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/"
    # base_dir+="SimOutput/ChainingCheckHigherModE"
    base_dir+="SimOutput/ChainingPhaseDiagramFinal"
    plot_dict = {}
    # run_nums=[256,259]

    if len(run_nums)!=1:
        quit()
    ax_lims=0

    min_min_x=np.inf
    max_max_x=-np.inf
    min_min_y=np.inf
    max_max_y=-np.inf
    for idx,ii in enumerate(run_nums):
        r=idx%2
        c=idx//2
        shift=np.array([r*500,c*350])
        file=glob(f"{base_dir}/run{ii}/repeat3/final*")[0]
        dat=pd.read_csv(file,sep='\t')
        dat['pos_x']+=shift[0]
        dat['pos_y']+=shift[1]

        maxx=dat['pos_x'].max()+2.5
        minx=dat['pos_x'].min()-2.5
        maxy=dat['pos_y'].max()+2.5
        miny=dat['pos_y'].min()-2.5
        max_max_x=max(max_max_x,maxx)
        min_min_x=min(min_min_x,minx)
        max_max_y=max(max_max_y,maxy)
        min_min_y=min(min_min_y,miny)
    ax_rng=max_max_x-min_min_x
    print(f"{ax_rng=}",f"{max_max_x=}",f"{min_min_x=}",f"{max_max_y=}",f"{min_min_y=}")
    files=[
        glob(f"{base_dir}/run{ii}/repeat3/final*.dat")[0]
        for ii in run_nums
        ]
    createImage(rn,files,bounds,ax_rng,scalebar=True)


    number_to_plot=max([
        len(os.listdir((f"{base_dir}/run{ii}/repeat3/"))) for ii in run_nums
        ])

    bounds=(min_min_x,max_max_x,min_min_y,max_max_y)
    for rn in reversed(range(number_to_plot)):
        try:
            files=[
                f"{base_dir}/run{ii}/repeat3/biofilm_{rn:05d}.dat"
                for ii in run_nums
                ]
            createImage(rn,files,bounds,ax_rng)
        except Exception as e:
            files=[
                glob(f"{base_dir}/run{ii}/repeat3/final*.dat")[0]
                for ii in run_nums
                ]
            createImage(rn,files,bounds,ax_rng)


    # for key,cels in plot_dict.items():
    #     base_dir="/home/rory/PhD_sync/BiofilmPhageDES/NBIC"
    #     save_name=f"{base_dir}/{key}_NBIC_plot.pdf"
    #     # lp=key*0.1
    #     # print(f"{lp=}")
    #     plotColony(cels,ax_lims,save_name,lp)

def plotColony(cells,ax_lim,save_name,lp):
    print(f"Creating {save_name}")

    fig,ax = plt.subplots(1,1,figsize=[20,20])
    ax.axis("scaled")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_xlim([-0.5*ax_lim,0.5*ax_lim])
    ax.set_ylim([-0.5*ax_lim,0.5*ax_lim])
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axis('off')
    # scalebar = AnchoredSizeBar(ax.transData,
    #                            100,
    #                            r'\SI{100}{\micro\meter}',
    #                            'lower right',
    #                            pad=0.1,
    #                            color='black',
    #                            frameon=False,
    #                            size_vertical=1,
    #                            fontproperties=fontprops)
    # ax.add_artist(scalebar)

    com=RodShapedBacterium.getCOM(cells)
    min_y=np.inf
    colour_fn=lambda th: cm.get_cmap('hsv')(getNematicAngle(th)/(np.pi))
    makeColourWheel(colour_fn)
    RodShapedBacterium.findMicrodomains(cells,colour_fn)
    for cell in cells:
        min_y=min(min_y,cell.rcm[1])
        cell.addElementToPlot(ax,colour=cell.colour)

    try:
        ax.text(com[0],max(-0.45*ax_lim,min_y-0.05*ax_lim),rf"$p_\ell={lp:1.1f}$",fontproperties=fontprops)
    except Exception as e:
        print(e)
        quit()
        ax.text(0,-0.45*ax_lim,rf"$p_\ell=$ {lp:1.1f}",fontproperties=fontprops)

    fig.savefig(save_name,
                format="pdf",
                bbox_inches='tight',
                dpi=dpi,
                transparent=True)

    plt.close()

if __name__ == "__main__":
    run_nums=[256]
    dpi=3000
    RodShapedBacterium.sus_vis_radius_factor=0.7
    createSingleImage(run_num=256,outnum=79)
    quit()
    makeComparisonMorphology()
