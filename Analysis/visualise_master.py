"""
    Script to provide helper functions to visualise phage from BiofilmDES

    This script sets up and creates a visualisation environment
"""

# Add parent to path
import sys

sys.path.insert(0, '..')

# Standard libraries
import matplotlib.pyplot as plt
from matplotlib.pyplot import axis, cm
from matplotlib import ticker
from matplotlib.collections import PatchCollection, LineCollection
import numpy as np
import numpy.linalg as linalg
import pandas as pd
import argparse
import os
import glob
import subprocess
import multiprocessing as mp
from joblib import Parallel, delayed
from copy import copy
import re
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
fontprops = fm.FontProperties(size=10)

# Third party
from shapely.geometry import Point
from hidden_prints import HiddenPrints
from descartes import PolygonPatch

# Custom modules
from DistributionFunctions import getColonyWavelength
from FindDefects import intialiseElementsFromData
from visualise_infected_biofilm import VisInfectedBiofilm
from visualise_biofilm import setFigSize,setClassStaticMembers
from phage_vis_class import Phage
from infected_vis_class import Infected
from RodShapedBacteria import RodShapedBacterium
from ChainingRodShapedBacteria import ChainingRodShapedBacterium as CRSB
from AG43RodShapedBacteria import AG43RodShapedBacterium
from SphericalBacteria import SphericalBacterium
from RodShapedBacteria import setProjection  # choose plane to visualise
from analysis_master import loadClusters,saveClusters
from fastCellPlotting import addAllCellsToPlot, addSpringLinks, addChainLinks
import utilities as ut

import matplotlib as mpl
ut.setMPL()

def getAxisLimits(file):
    dat=pd.read_csv(file,sep='\t')
    maxx=dat['pos_x'].max()+2.5
    minx=dat['pos_x'].min()-2.5
    maxy=dat['pos_y'].max()+2.5
    miny=dat['pos_y'].min()-2.5
    ax_rng=max(maxx-minx,maxy-miny)
    return ax_rng

def getMaxRng(file):
    dir=os.path.split(file)[0]
    use_file=glob.glob(f"{dir}/final*")[0]
    return getAxisLimits(use_file)

class AnyObject(object):
    pass

def addColourWheel(fig,colour_fn):

    # print("load colour wheel from file to add")
    # quit()
    left, bottom, width, height = [0.7, 0.825, 0.15, 0.15]
    ax = fig.add_axes([left, bottom, width, height],polar=True)

    angles = np.linspace(-np.pi, np.pi, 360)
    for sa,ea in zip(angles[:-1],angles[1:]):
        ax.fill_between([sa,ea],0,1,color=colour_fn(sa))
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_axis_off()

def getNematicAngle(alpha):
    if np.abs(alpha)>np.pi:
        print(f"alpha is too big! {alpha:2.3f}")
        quit()

    if alpha>=0:
        return alpha

    if alpha<0:
        return np.pi+alpha

def makeMovie(start,
              end,
              step=1,
              max_cores=mp.cpu_count(),
              base_path="../GeneratedOutput",
              data_dir="SimOutput",
              animation_dir="animation",
              stem="biofilm",
              class_dict={
                 "RodShaped"         : RodShapedBacterium,
                 "Spherical"         : RodShapedBacterium,
                 "ChainingRodShaped" : CRSB
                 },
              track=False):
    """
    Super helpful description
    """
    # print("Make a better movie function")
    # quit()
    # Convert to parallel
    # for ii in range(args.start_index,args.end_index+1):
    #     args.sus_file = f"data/biofilm_{ii:05d}.dat"
    #
    #     args.inf_file = f"data/main_infected_{ii:05d}.dat"
    #
    #     args.phg_file = f"data/main_phage_{ii:05d}.dat"
    #
    #     args.fig_file = f"output/infected_biofilm_{ii:05d}.png"
    #
    #     createFull2DPlot(
    #         sus_file=args.sus_file,inf_file=args.inf_file,
    #         phg_file=args.phg_file,fig_file=args.fig_file,
    #         colour=cm.viridis,show=False
    #         )
    try:
        log_file = f"{base_path}/{data_dir}/{stem}.log"
        hps=pd.read_csv(log_file,sep="\t").to_dict('list')
        bending_moduli=hps['BendRig'][0]
        linking_prob=hps['LinkingProb'][0]
    except Exception as e:
        print(e)
        log_file = f"{base_path}/{data_dir}/World.log"
        hps=pd.read_csv(log_file,sep="\t").to_dict('list')
        # bending_moduli=hps['BendRig'][0]
        # linking_prob=hps['LinkingProb'][0]
        print(hps)
        bending_moduli=0
        linking_prob=0

    sus_data = pd.read_csv(f"{base_path}/{data_dir}/{stem}_{0:05d}.dat", sep='\t')
    all_sus, species_load_dict = intialiseElementsFromData(sus_data,class_dict)
    ax_lim = np.max(VisInfectedBiofilm.getAxisLimits(all_sus)) + 5
    print(ax_lim)

    if track==True:
        ax_lim_max=getMaxRng(f"{base_path}/{data_dir}/{stem}_{0:05d}.dat")
    else:
        ax_lim_max=None

    """
        Generate list of tuples for each input to createFull2DPlot in the
        order
        sus_file,inf_file,phg_file,fig_file,sus_class=RodShapedBacterium,
        inf_class=Infected,phg_class=Phage,plane_loc=None,eps=1e-4,
        projection="xy",colour=cm.Greens,show=False
    """
    if end==-1:
        end=len(glob.glob(f"{base_path}/{data_dir}/{stem}_*.dat"))-1

    print(f"Running with {max_cores} cores")
    par_args_list = [
        {
            "in_file": f"{base_path}/{data_dir}/{stem}_{ii:05d}.dat",
            "fig_file": f"{base_path}/{animation_dir}/{stem}_{ii:05d}.png",
            "colour": cm.viridis,
            "show": False,
            "ax_lim": ax_lim_max,
            "dpi": 600,
            "class_dict" : class_dict
        } for ii in reversed(range(start, end + 1, step))]
    Parallel(n_jobs=max_cores, verbose=10, backend="multiprocessing")(
        delayed(createFull2DPlot)(**par_args)
        for par_args in par_args_list
        )
    movie_command = (
            f"ffmpeg -y -framerate 10 -i {base_path}/{animation_dir}/{stem}_%05d.png"
            + " -c:v copy -vb 20M"
            + f" {base_path}/{animation_dir}/{stem}_{linking_prob=}_{bending_moduli=}.mp4"
            )
    os.system(movie_command)
    quit()

    dr=0.1
    par_args_list = [
        {
            "in_file": f"{base_path}/{data_dir}/{stem}_{ii:05d}.dat",
            "fig_file": f"{base_path}/{animation_dir}/Q/{stem}_Q_{jj:05d}.png",
            "dr" : dr,
            "colour": cm.viridis,
            "show": False,
            "ax_lim": ax_lim,
            "dpi": 200,
            "class_dict" : class_dict
        } for jj,ii in enumerate(reversed(range(start, end + 1, 10)))]
    Parallel(n_jobs=max_cores, verbose=10, backend="multiprocessing")(
        delayed(createOrderParameterPlot)(**par_args)
        for par_args in par_args_list
        )
    movie_command = (
            f"ffmpeg -y -framerate 5 -i {base_path}/{animation_dir}/Q/{stem}_defect_{dr=}_%05d.png"
            + " -c:v copy -vb 20M"
            + f" {base_path}/{animation_dir}/Q/{stem}_defect.mp4"
            )
    os.system(movie_command)
    quit()

    # Convert to Parallel

    # Parallel(n_jobs=max_cores, verbose=10, backend="multiprocessing")(
    #     delayed(
    #         suppressOutput(createFull2DPlot))(**par_args)
    #     for par_args in par_args_list
    # )

    par_args_list = [
        {
            "in_file": f"{base_path}/{data_dir}/{stem}_{ii:05d}.dat",
            "fig_file": f"{base_path}/{animation_dir}/density/{stem}_density_{ii:05d}.png",
            "colour": cm.viridis,
            "show": False,
            "ax_lim": ax_lim,
            "dpi": 300,
            "class_dict" : class_dict
        } for ii in reversed(range(start, end + 1, step))]
    Parallel(n_jobs=max_cores, verbose=10, backend="multiprocessing")(
        delayed(createDensityPlot)(**par_args)
        for par_args in par_args_list
        )
    movie_command = (
            f"ffmpeg -y -framerate 10 -i {base_path}/{animation_dir}/density/{stem}_density_%05d.png"
            + " -c:v copy -vb 20M"
            + f" {base_path}/{animation_dir}/density/{stem}_density.mp4"
            )
    os.system(movie_command)


    par_args_list = [
        {
            "in_file": f"{base_path}/{data_dir}/{stem}_{ii:05d}.dat",
            "fig_file": f"{base_path}/{animation_dir}/energy/{stem}_energy_{ii:05d}.png",
            "colour": cm.viridis,
            "show": False,
            "ax_lim": ax_lim,
            "dpi": 200,
            "class_dict" : class_dict
        } for ii in reversed(range(start, end + 1, step))]
    Parallel(n_jobs=max_cores, verbose=10, backend="multiprocessing")(
        delayed(createEnergyPlot)(**par_args,energy_type='both')
        for par_args in par_args_list
        )


    # with mp.Pool(processes=args.max_cores) as pool:
    #     results = pool.starmap(
    #         createFull2DPlot,
    #         [par_args for par_args in par_args_list]
    #         )

    # if step == 1:

    # else:
    #     movie_command = (
    #             "ffmpeg -y -f image2 -pattern_type glob -r 10"
    #             + f" -i '{base_path}/{animation_dir}/{stem}_*.png'"
    #             + " -c:v copy -vf fps=25 -vb 20M"
    #             + f" {base_path}/{animation_dir}/{stem}.mp4"
    #     )
    # os.system(movie_command)
    # os.system(f"ffplay {base_path}/{animation_dir}/{stem}.mp4")
    movie_command = (
            f"ffmpeg -y -framerate 10 -i {base_path}/{animation_dir}/energy/{stem}_energies_%05d.png"
            + " -c:v copy -vb 20M"
            + f" {base_path}/{animation_dir}/energy/{stem}_both.mp4"
            )
    os.system(movie_command)
    # os.system(f"ffplay {base_path}/{animation_dir}/energy/{stem}.mp4")

def suppressOutput(func):
    def noPrint(*args, **kwargs):
        with HiddenPrints():
            out = func(*args, **kwargs)
        return out

    return noPrint

def createFull2DPlot(in_file,
                     fig_file,
                     class_dict={
                        "RodShaped"         : RodShapedBacterium,
                        "Spherical"         : RodShapedBacterium,
                        "ChainingRodShaped" : CRSB
                        },
                     plane_loc=None,
                     eps=1e-4,
                     projection="xy",
                     colour=cm.Greens,
                     show=False,
                     colour_chains=False,
                     dpi=200,
                     ax_lim=None,
                     transparent=False,
                     include_time=True):
    """
        Creates a plot of a slice through the biofilm centered on the chosen
        location, within tolerance (eps) in the given plane (projection)

        Parameters:
            in_file: str
                Data file expected with column headers
                cell_type,cell_id,length,radius,pos_x,pos_y,pos_z,ori_x,ori_y,ori_z
                with '\t' delimeters
            class_dict: dict
                These define how objects will appear on the axes. The list should
                be a dictionary of cell_type keys with value the class to use to
                represent this.
            fig_file: str
                This is the fig_file  image is saved to.
            plane_loc: float
                value in the projected direction to center the plots on
                e.g. if projection="xy" then this is the z direction
            eps: float
                bacteria with centers of mass in |plane_loc - com[select]|<eps
                will be displayed where select picks the direction of projection
                e.g. if projection="xy" then select = 2 for the z direction
            projection: str
                determines onto which plane to project the 3D biofilm
                options: "xy","xz","yz"
            colour: callable -> RGBA
                colourmap function
            show: bool
                If True show the figure after saving and any outputs to terminal
            dpi: int
                dots per inch - 300 is high quality but takes a long time
            ax_lim: float
                axis range with x,y in [-0.5 ax_lim, 0.5 ax_lim]
            include_time: bool
                if true, try to add the time to the plot corresponds to

        Effect:
            Saves the slice
    """
    add_links=False # Add chain links to plot
    bk='w'
    if bk=='w':
        fg='k'
    elif bk=='k':
        fg='w'
    # If this is True, then the background will not be black
    transparent=False
    # Set up axes and load susceptible
    if os.path.exists(fig_file):
        print(f"{fig_file} already exists, skipping")
        return
    try:
        vis_bio = VisInfectedBiofilm(in_file,class_dict)
    except Exception as e:
        raise
    vis_bio.setSaveName(fig_file)
    vis_bio.setUpFigAx(fig_size=setFigSize(2*247),bk=None)
    # If left blank, use class default colour schemes
    # vis_bio.setElementColourScheme()

    # ------------------------ Quasi 2D considerations -------------------------
    tangent_axes, normal_axes = setProjection(projection)

    # Select cells to display which have an interesction wihthin eps of the
    # plane location
    if plane_loc is not None:
        vis_bio.element_list = [
            element for element in vis_bio.element_list
            if cellPlaneIntersect(element, plane_loc, eps, normal_axes)
        ]

        # Sort cells so image appears as it would along the projected line of sight
        vis_bio.element_list = sorted(vis_bio.element_list,
                                      key=lambda cell: cell.rcm[normal_axes]
                                      )
    # ---------------------------------------------------------------------------

    # Choose axis size so that all elements are displayed

    if ax_lim==None:
        ax_lim=getAxisLimits(in_file)
    # else:
    #     ax_rng=getAxisLimits(in_file)
    #     ax_rng=2**np.ceil(np.log2(ax_rng))
    #     ax_rng+=5
    #     ax_lim=min(ax_lim,ax_rng)
    # The ax_rng is required to make the edge of the cells the right thickness

    colour_wheel=True
    if colour_wheel:
        addColourWheel(vis_bio.fig,RodShapedBacterium.colour_fn)

    cluster_flag=False
    if cluster_flag:
        cluster_file=vis_bio.save_name.replace("figures","AnalysisResults")
        cluster_file=cluster_file.replace(".pdf","_cluster.dat")
        try:
            clusters=loadClusters(vis_bio.element_list,cluster_file)
            for key,cluster in clusters.items():
                cc=RodShapedBacterium.colour_fn(cluster[0].theta)
                for cell in cluster:
                    cell.colour=cc
        except Exception as e:
            print(e)
            clusters=RodShapedBacterium.findMicrodomains(
                vis_bio.element_list,
                RodShapedBacterium.colour_fn
                )
            saveClusters(clusters,cluster_file)
        overlaps=RodShapedBacterium.getOverlaps(vis_bio.element_list)
        print(f"maximum {np.max(overlaps)=}")
        overlap_file=vis_bio.save_name.replace("figures","AnalysisResults")
        overlap_file=overlap_file.replace(".pdf","_overlap.pdf")
        fig,ax=plt.subplots()
        ax.hist(overlaps,density=True)
        ax.set_xlabel("overlap")
        fig.savefig(overlap_file)
    else:
        for cell in vis_bio.element_list:
            cell.colour=RodShapedBacterium.colour_fn(cell.theta)

    # try:
    #     AG43RodShapedBacterium.createSpringHashList(vis_bio.element_list)
    # except Exception as e:
    #     print(e)
    #     quit()
    # for cell in vis_bio.element_list:
    #     cell.addElementToPlot(vis_bio.ax, ax_rng=ax_lim, show_id=False,
    #                           colour=cell.colour
    #                           )

    RodShapedBacterium.makeHashList(vis_bio.element_list)
    # if colour_chains:
    #     print("colouring chains")
    #     chains=CRSB.getChains(vis_bio.element_list)
        # norm = plt.Normalize(0,0.5*np.pi)
        # for chain in chains:
        #     thetas=CRSB.getCurvature(chain)
        #     # if len(thetas):
        #     #     thetas=np.convolve(thetas,np.full(2,0.5),mode='full')
        #     #     thetas[0]*=2
        #     #     thetas[-1]*=2
        #     # else:
        #     #     thetas=[0]
        #     try:
        #         chain_colour=cm.get_cmap('PuOr')(np.max(thetas))
        #     except Exception as e:
        #          chain_colour=cm.get_cmap('PuOr')(0)
        #     for ii,cell in enumerate(chain):
        #         # cell.colour=cm.get_cmap('viridis')(thetas[ii])
        #         cell.colour=chain_colour

    # wl,th,rcms=getColonyWavelength(vis_bio.element_list,nbins=100)
    addAllCellsToPlot(
        cells=vis_bio.element_list,
        ax=vis_bio.ax,
        ax_rng=ax_lim,
        show_id=False,
        ec=fg
        )
    if add_links:
        addChainLinks(vis_bio.element_list,ax=vis_bio.ax,
                      ax_rng=ax_lim,show_id=False,colour=fg
                      )


    # f=lambda x,a,x0,lam: a*np.sin(2*np.pi*(x-x0)/lam)
    # vis_bio.ax.scatter(rcms[:,0],rcms[:,1],s=20,fc='w',ec='k')
    # vis_bio.ax.plot(rcms[:,0],f(rcms[:,0],*wl))

    # Make this one function
    # try:
    #     addChainLinks(
    #         cells=vis_bio.element_list,
    #         ax=vis_bio.ax,
    #         ax_rng=ax_lim,
    #         show_id=False
    #         )
    # except Exception as e:
    #     try:
    #         addSpringLinks(
    #             cells=vis_bio.element_list,
    #             ax=vis_bio.ax,
    #             ax_rng=ax_lim,
    #             show_id=False
    #             )
    #     except Exception as e:
    #         print(e)

    # Create a nice legend
    legend_flag=False
    include_time=False
    if legend_flag:
        legend_labels=[]
        legend_elements=[]
        handler_map={}
        for key,val in vis_bio.species_load_dict.items():
            legend_labels.append(f"{key} : {val}")
            Class=class_dict[key]
            try:
                try:
                    cell=Class(0,0,0,0,0,0,1,0,0)
                except Exception as e:
                    cell=Class(0,0,0,0,0,0,1,0,0,'')
            except Exception as e:
                cell=Class(0,0,0,0,0)
            cell.colour=RodShapedBacterium.colour_fn(0)
            legend_elements.append(cell)
            handler_map[Class]=cell
        if include_time:
            try:
                print("Read in log file for output frequency")
                quit()
                time=int(re.findall(r'\d+',os.path.split(in_file)[1])[0])*0.02
                title=f"Time: {time:2.2f} hrs"
            except Exception as e:
                print(e)
                time=None
                title=None
        else:
            time=None
            title=None

        vis_bio.ax.legend(legend_elements,
                          legend_labels,
                          handler_map=handler_map,
                          title=title,
                          loc='upper right', bbox_to_anchor=(1.1,1.1),
                          fontsize='small',
                          title_fontsize='small',
                          borderaxespad=0,
                          labelspacing=0.5,
                          handlelength=2.5,
                          borderpad=1.25,
                          framealpha=0,
                          frameon=False)

    # Create a scale bar
    sc=10**int(np.floor(np.log10(ax_lim)))
    scalebar = AnchoredSizeBar(vis_bio.ax.transData,
                               sc,
                               '\SI{{{:d}}}'.format(sc)+'{\micro\meter}',
                               'lower right',
                               pad=0.1,
                               color=fg,
                               frameon=False,
                               size_vertical=0,
                               alpha=0.7,
                               fontproperties=fontprops)
    vis_bio.ax.add_artist(scalebar)

    # Need this so that the polygons look nice
    vis_bio.ax.axis("scaled")
    vis_bio.ax.axes.get_xaxis().set_visible(False)
    vis_bio.ax.axes.get_yaxis().set_visible(False)
    vis_bio.ax.axis('off')

    # projection in xy is approximately circular
    vis_bio.ax.set_xlim([-0.5 * ax_lim, 0.5 * ax_lim])
    vis_bio.ax.set_ylim([-0.5 * ax_lim, 0.5 * ax_lim])

    if projection != "xy":
        vis_bio.ax.set_ylim([None, None])

    vis_bio.fig.savefig(vis_bio.save_name,
                        format=os.path.splitext(fig_file)[1].replace('.', ''),
                        bbox_inches='tight',
                        transparent=transparent,
                        dpi=dpi)
    if show == True:
        plt.show()
    plt.close()

def createEnergyPlot(
    in_file,
    fig_file,
    energy_type,
    class_dict={
        "RodShaped"         : RodShapedBacterium,
        "Spherical"         : RodShapedBacterium,
        "ChainingRodShaped" : CRSB
    },
    plane_loc=None,
    eps=1e-4,
    projection="xy",
    colour=cm.Greens,
    show=False,
    dpi=200,
    bend_only=False,
    ax_lim=None,
    transparent=False,
    include_time=True):

    if energy_type=='Hertzian':
        fig_file=fig_file.replace('_energy_','_contact_energy_')
    elif energy_type=='Chaining':
        fig_file=fig_file.replace('_energy_',f'_chaining_energy_{bend_only=}_')
    elif energy_type=='both':
        fig_file=fig_file.replace('_energy_',f'_energies_{bend_only=}_')
    else:
        print(f"{energy_type=}")
        quit()
    print(f"{fig_file=}")

    if os.path.exists(fig_file):
        print(f"{fig_file} already exists, skipping")
        return

    try:
        element_list, species_load_dict = intialiseElementsFromData(
            pd.read_csv(in_file,sep='\t'),
            class_dict
            )
    except Exception as e:
        print(e)
        raise

    RodShapedBacterium.getContactEnergies(element_list)

    if energy_type=='Hertzian':
        energy_field, energy_grid = RodShapedBacterium.computeEnergyGrid(
            element_list,
            dr=0.5
            )
        # fig_file=fig_file.replace('_energy_','_contact_energy_')
    elif energy_type=='Chaining':
        energy_field, energy_grid = CRSB.computeEnergyGrid(
            element_list,
            dr=0.5,
            bend_only=bend_only
            )
        # fig_file=fig_file.replace('_energy_','_chaining_energy_')
    elif energy_type=='both':
        energy_field, energy_grid = RodShapedBacterium.computeEnergyGrid(
            element_list,
            dr=0.5
            )
        spring_field, spring_grid = CRSB.computeEnergyGrid(
            element_list,
            dr=0.5,
            bend_only=bend_only
            )
        # fig_file=fig_file.replace('_energy_','_energies_')
    else:
        print(f"{energy_type=}")
        quit()
    print(f"{fig_file=}")

    if bend_only:
        spring_label="Bend energy"
    else:
        spring_label="Total spring energy"

    energy_field[energy_field<1e-12]=np.nan
    if energy_type=='both':
        spring_field[spring_field<1e-12]=np.nan
        fig,axs=plt.subplots(2,1,figsize=[20,10])
        for ax in axs:
            ax.axis("scaled")
        ax=axs[0]
    else:
        fig,ax=plt.subplots(1,1,figsize=[10,10])
        ax.axis("scaled")

    lev_exp = np.arange(-3,1,0.5,dtype=float)
    levs = np.power(10, lev_exp)

    # levels=np.arange(np.min(energy_field),np.max(energy_field),endpoint=True)
    if energy_type!='Chaining':
        CS=ax.contourf(*energy_grid,energy_field,alpha=0.8,
                       cmap='inferno',zorder=2,
                       locator=ticker.LogLocator())
    else:
        CS=ax.contourf(*energy_grid,energy_field,alpha=0.8,
                       # levels=levs,norm=mpl.colors.LogNorm(),
                       cmap='inferno',zorder=2,
                       locator=ticker.LogLocator())
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)


    X,Y = energy_grid
    # ax.set_xlim([np.min(X),np.max(X)])
    # ax.set_ylim([np.min(Y),np.max(Y)])
    cbar=fig.colorbar(CS,cax=cax)
    if energy_type!='Chaining':
        cbar.set_label("Contact energy",size=14)
        cbar.ax.tick_params(labelsize=12)
    else:
        cbar.set_label(spring_label,size=14)
        cbar.ax.tick_params(labelsize=12)


    for cell in element_list:
        cell.colour="None"

    if energy_type=='both':
        ax=axs[1]
        # levels=np.arange(np.min(energy_field),np.max(energy_field),endpoint=True)
        CS=ax.contourf(*spring_grid,spring_field,alpha=0.8,
                       # levels=levs,norm=mpl.colors.LogNorm(),
                       cmap='inferno',zorder=2,
                       locator=ticker.LogLocator())
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        X,Y = spring_grid
        # ax.set_xlim([np.min(X),np.max(X)])
        # ax.set_ylim([np.min(Y),np.max(Y)])
        cbar=fig.colorbar(CS,cax=cax)
        cbar.set_label(spring_label,size=14)
        cbar.ax.tick_params(labelsize=12)

    filedir = os.path.dirname(fig_file)
    # print(filedir)
    if filedir == '':
        pass
    elif not os.path.exists(filedir):
        os.makedirs(filedir)

    try:
        for ax in axs:
            ax.set_yticks([])
            ax.set_xticks([])
    except Exception as e:
        ax.set_yticks([])
        ax.set_xticks([])

    # Choose axis size so that all elements are displayed
    if ax_lim==None:
        ax_lim = 0
        try:
            ax_lim = max(VisInfectedBiofilm.getAxisLimits(element_list), ax_lim)
            ax_lim = 2**np.ceil(np.log2(ax_lim))
        except Exception as e:
            print(e)
            quit()
    else:
        ax_lim=max(np.max(X)-np.min(X),np.max(Y)-np.min(Y))

    try:
        for ax in axs:
            ax.set_xlim([-0.5 * ax_lim, 0.5 * ax_lim])
            ax.set_ylim([-0.5 * ax_lim, 0.5 * ax_lim])
            addAllCellsToPlot(
                cells=element_list,
                ax=ax,
                ax_rng=ax_lim,
                alpha=0.5,
                show_id=False
                )
    except Exception as e:
        ax.set_xlim([-0.5 * ax_lim, 0.5 * ax_lim])
        ax.set_ylim([-0.5 * ax_lim, 0.5 * ax_lim])
        addAllCellsToPlot(
            cells=element_list,
            ax=ax,
            ax_rng=ax_lim,
            alpha=0.5,
            show_id=False
            )

    fig.savefig(fig_file,
                format=os.path.splitext(fig_file)[1].replace('.', ''),
                bbox_inches='tight',
                transparent=transparent,
                dpi=dpi)

def createDensityPlot(
    in_file,
    fig_file,
    class_dict={
        "RodShaped"         : RodShapedBacterium,
        "Spherical"         : RodShapedBacterium,
        "ChainingRodShaped" : CRSB
    },
    plane_loc=None,
    eps=1e-4,
    projection="xy",
    colour=cm.Greens,
    show=False,
    dpi=200,
    ax_lim=None,
    transparent=False,
    include_time=True):

    print(f"{fig_file=}")

    if os.path.exists(fig_file):
        print(f"{fig_file} already exists, skipping")
        return

    try:
        element_list, species_load_dict = intialiseElementsFromData(
            pd.read_csv(in_file,sep='\t'),
            class_dict
            )
    except Exception as e:
        print(e)
        raise

    density_field, density_grid = RodShapedBacterium.computeDensity(
        element_list,
        dr=0.5
        )

    density_field[density_field<1e-5]=np.nan
    fig,ax=plt.subplots(1,1,figsize=[10,10])
    ax.axis("scaled")

    CS=ax.contourf(*density_grid,density_field,alpha=0.8,
                   # levels=levs,norm=mpl.colors.LogNorm(),
                   cmap='inferno',zorder=2,
                   locator=ticker.LogLocator())
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)


    X,Y = density_grid
    # ax.set_xlim([np.min(X),np.max(X)])
    # ax.set_ylim([np.min(Y),np.max(Y)])
    cbar=fig.colorbar(CS,cax=cax)
    cbar.set_label(r"$\rho$",size=14)
    cbar.ax.tick_params(labelsize=12)


    for cell in element_list:
        cell.colour="None"

    filedir = os.path.dirname(fig_file)
    # print(filedir)
    if filedir == '':
        pass
    elif not os.path.exists(filedir):
        os.makedirs(filedir)

    ax.set_yticks([])
    ax.set_xticks([])

    # Choose axis size so that all elements are displayed
    if ax_lim==None:
        ax_lim = 0
        try:
            ax_lim = max(VisInfectedBiofilm.getAxisLimits(element_list), ax_lim)
            ax_lim = 2**np.ceil(np.log2(ax_lim))
        except Exception as e:
            print(e)
            quit()
    else:
        ax_lim=max(np.max(X)-np.min(X),np.max(Y)-np.min(Y))

    ax.set_xlim([-0.5 * ax_lim, 0.5 * ax_lim])
    ax.set_ylim([-0.5 * ax_lim, 0.5 * ax_lim])
    addAllCellsToPlot(
        cells=element_list,
        ax=ax,
        ax_rng=ax_lim,
        alpha=0.5,
        show_id=False
        )

    fig.savefig(fig_file,
                format=os.path.splitext(fig_file)[1].replace('.', ''),
                bbox_inches='tight',
                transparent=transparent,
                dpi=dpi)

def createOrderParameterPlot(
    in_file,
    fig_file,
    dr=0.1,
    class_dict={
        "RodShaped"         : RodShapedBacterium,
        "Spherical"         : RodShapedBacterium,
        "ChainingRodShaped" : CRSB
    },
    plane_loc=None,
    eps=1e-4,
    projection="xy",
    colour=cm.Greens,
    show=False,
    dpi=200,
    ax_lim=None,
    transparent=False,
    include_time=True):

    print(f"{fig_file=}")
    fig_file=fig_file.replace('_Q_',f'_Q_{dr=}_')
    datadir = os.path.dirname(fig_file).replace("figures","data")
    datadir = datadir.replace("animation","data")
    print(os.path.splitext(fig_file))
    ext = os.path.splitext(in_file)[-1]
    fname = os.path.basename(in_file).replace(ext,'.npy').replace('biofilm_','')
    data_file=f"{datadir}/Q_{dr=}_{fname}"
    if os.path.exists(fig_file):
        print(f"{fig_file} already exists, skipping")
        return
    if not os.path.exists(datadir):
        os.makedirs(datadir)

    try:
        element_list, species_load_dict = intialiseElementsFromData(
            pd.read_csv(in_file,sep='\t'),
            class_dict
            )
    except Exception as e:
        print(e)
        raise

    # print("Change to work on chunks insteas (check np splits)")
    # print("Also check if what the director field looks like")
    # quit()
    q_field, Q_field, Q_grid = RodShapedBacterium.computeChargeDensity(
        element_list,
        dr=dr,
        fname=data_file
        )

    # for ii,l in enumerate(locs):
    #     q_min = q_field[r]
    #     r_min = X[r],Y[r]
    #     print(f"{q_min} @ {r_min}")

    # dr=2
    # q_field, Q_field, Q_grid = RodShapedBacterium.computeChargeDensityWinding(
    #     element_list,
    #     dr=dr
    #     )

    fig,ax=plt.subplots(1,1,figsize=[10,10])
    ax.axis("scaled")
    cmap='RdYlGn'

    q_plot_field=np.copy(q_field)
    q_plot_field[np.abs(q_plot_field)<1e-5]=np.nan
    CS=ax.contourf(*Q_grid,q_plot_field,alpha=0.8,
                   # levels=levs,norm=mpl.colors.LogNorm(),
                   cmap=cmap,zorder=2)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    X,Y = Q_grid

    # ax.set_xlim([np.min(X),np.max(X)])
    # ax.set_ylim([np.min(Y),np.max(Y)])
    cbar=fig.colorbar(CS,cax=cax)
    cbar.set_label(rf"$q$",size=14)
    cbar.ax.tick_params(labelsize=12)

    # for cell in element_list:
    #     cell.colour="None"

    filedir = os.path.dirname(fig_file)
    # print(filedir)
    if filedir == '':
        pass
    elif not os.path.exists(filedir):
        os.makedirs(filedir)

    ax.set_yticks([])
    ax.set_xticks([])

    # Choose axis size so that all elements are displayed
    if ax_lim==None:
        ax_lim = 0
        try:
            ax_lim = max(VisInfectedBiofilm.getAxisLimits(element_list), ax_lim)
            ax_lim = 2**np.ceil(np.log2(ax_lim))
        except Exception as e:
            print(e)
            quit()
    else:
        ax_lim=max(np.max(X)-np.min(X),np.max(Y)-np.min(Y))

    ax.set_xlim([-0.5 * ax_lim, 0.5 * ax_lim])
    ax.set_ylim([-0.5 * ax_lim, 0.5 * ax_lim])
    addAllCellsToPlot(
        cells=element_list,
        ax=ax,
        ax_rng=0.5*ax_lim,
        alpha=0.3,
        show_id=False
        )

    fig.savefig(fig_file,
                format=os.path.splitext(fig_file)[1].replace('.', ''),
                bbox_inches='tight',
                transparent=transparent,
                dpi=dpi)


    fig,ax=plt.subplots(1,1,figsize=[10,10])
    ax.axis("scaled")

    ax.set_yticks([])
    ax.set_xticks([])

    ax.set_xlim([-0.5 * ax_lim, 0.5 * ax_lim])
    ax.set_ylim([-0.5 * ax_lim, 0.5 * ax_lim])
    addAllCellsToPlot(
        cells=element_list,
        ax=ax,
        ax_rng=0.5*ax_lim,
        alpha=0.3,
        show_id=False
        )

    try:
        # cfn = cm.get_cmap(cmap)
        X,Y = Q_grid

        X_mins,X_maxs=RodShapedBacterium.filterDefects(q_field,Q_grid)
        charge,pos=RodShapedBacterium.findLocalDefects(
            np.vstack((X_mins,X_mins)),
            element_list
            )
        # colour = cfn(0.5*(1+np.sign(np.min(q_mins))))
        sc=ax.scatter(
            pos[:,0],pos[:,1],
            c=charge,
            cmap=cmap,
            zorder=4,alpha=0.7
            )
    except Exception as e:
        print(e)
        norm = mpl.colors.Normalize(vmin=-0.5,vmax=0.5)
        sc = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sc.set_array([])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar=fig.colorbar(sc,cax=cax)
    cbar.set_label(rf"$q$",size=14)
    cbar.ax.tick_params(labelsize=12)

    fig.savefig(fig_file.replace('_Q_',f'_defect_'),
                format=os.path.splitext(fig_file)[1].replace('.', ''),
                bbox_inches='tight',
                transparent=transparent,
                dpi=dpi)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='visualise_phage.py',
        usage='%(prog)s [options] path',
        description='Visualisation script for infected biofilms',
        epilog='Author: Rory Claydon'
    )
    parser.add_argument("-S", "--start-index", dest="start_index",
                        type=int,
                        required=True,
                        help="load data from"
                             + " SimOutput/output_{start_index:05d}.dat")
    parser.add_argument("-E", "--end-index",dest="end_index",
                        type=int,
                        required=False,
                        default=None,
                        help="load data up until"
                             + " SimOutput/output_{start_index:05d}.dat"
                             + " if end=-1 use all in directory")
    parser.add_argument("-P","--phage",action="store_true",required=False)
    parser.add_argument("--step", type=int, required=False,
                        default=1,
                        help="skip step number outputs")
    parser.add_argument("-I", "--in-file", dest="in_file",
                        type=str,
                        required=False,
                        default=None,
                        help="output data filename if different from default")
    parser.add_argument("-F", "--fig-file",dest="fig_file", type=str, default=None,
                        required=False, help="output figure filename")
    parser.add_argument("--max-cores", type=int, default=mp.cpu_count(),
                        required=False, action="store", dest="max_cores",
                        help="maximum number of cores to use to produce output")
    parser.add_argument("--dpi", type=int, default=200,
                        required=False, action="store", dest="dpi",
                        help="png resolution in dpi")
    parser.add_argument("--low-res", dest="low_res", action="store_true", default=False,
                        help="use lower resolution for faster png generation")
    parser.add_argument("--high-res", dest="high_res", action="store_true", default=False,
                        help="use higher resolution for better movie generation")
    parser.add_argument("-BD", "--base-dir", type=str, required=False, default="../GeneratedOutput",
                        action="store", dest="base_path",
                        help="default directory to look for inputs")
    parser.add_argument("-DD", "--data-dir", type=str, required=False, default="SimOutput",
                        action="store", dest="data_dir",
                        help="default directory name in base for inputs")
    parser.add_argument("-AD", "--anim-dir", type=str, required=False, default="animation",
                        action="store", dest="animation_dir",
                        help="default directory name in base for inputs")
    parser.add_argument("--stem", type=str, required=False, default="biofilm",
                        action="store", dest="stem",
                        help="default stem for input file names i.e. path/{stem}_%05d.dat")
    parser.add_argument("--log_file", type=str, required=False, default="biofilm",
                        action="store", dest="log_file.log",
                        help="default stem for input file names i.e. path/{log_file}")
    parser.add_argument("--fig-format", type=str, required=False, default="png",
                        action="store", dest="fig_format",
                        help="set figure format")
    parser.add_argument("--sus-vis-radius-factor", type=float, required=False, default=1,
                        action="store", dest="sus_vis_radius_factor",
                        help="Set visualisation width this fraction of real width for bacteria")
    parser.add_argument("--colour-chains",required=False,
                        action="store_true", dest="colour_chains",
                        help="colour chains instead of by angle")
    args = parser.parse_args()

    RodShapedBacterium.sus_vis_radius_factor=args.sus_vis_radius_factor

    try:
        log_file = f"{args.base_path}/{args.data_dir}/{args.log_file}"
        hps=setClassStaticMembers(log_file)
    except Exception as e:
        print(e)
        print("Continuing without log file information")

    if args.high_res == True and args.low_res == True:
        print("Please select one option of --low-res or --high-res")
        print(f"Using default dpi={args.dpi}")
    if args.high_res == False and args.low_res == False:
        print(f"Using dpi={args.dpi}")
    if args.low_res:
        args.dpi = 100
        print(f"Creating low resolution plots; dpi={args.dpi}")
    if args.high_res:
        args.dpi = 300
        print(f"Creating high resolution plots; dpi={args.dpi}")

    class_dict={
       "RodShaped"         : RodShapedBacterium,
       "Spherical"         : RodShapedBacterium,
       "ChainingRodShaped" : CRSB,
       "AG43RodShaped"     : AG43RodShapedBacterium
       }
    if args.phage:
        class_dict["Phage"]=Phage

    if args.start_index==-1:
        start_index=len(glob.glob(f"{args.base_path}/{args.data_dir}/{args.stem}_*.dat"))-1
    else:
        start_index=args.start_index

    # this logic is horrible
    if args.end_index == None:
        if args.in_file == None:
            # if there is only one file and we want the last output use final
            if args.start_index==-1:
                start_index+=1
                args.in_file = f"{args.base_path}/{args.data_dir}/final_{start_index:05d}.dat"
                args.stem="final"
            else:
                args.in_file = f"{args.base_path}/{args.data_dir}/{args.stem}_{start_index:05d}.dat"

        if args.fig_file == None:
            args.fig_file = f"{args.base_path}/{args.animation_dir}/{args.stem}_{start_index:05d}.{args.fig_format}"
            args.energy_fig_file = f"{args.base_path}/{args.animation_dir}/energy/{args.stem}_energy_{start_index:05d}.{args.fig_format}"


        print(f"Reading from output file {args.in_file}")
        print(f"Writing to figure file {args.fig_file}")
        createFull2DPlot(
            in_file=args.in_file,
            fig_file=args.fig_file,
            colour=cm.viridis,
            class_dict=class_dict,
            show=False,
            ax_lim=None,
            colour_chains=args.colour_chains,
            include_time=False
            )
        quit()
        try:
            createFull2DPlot(
                in_file=args.in_file,
                fig_file=args.fig_file,
                colour=cm.viridis,
                class_dict=class_dict,
                show=False,
                ax_lim='full',
                include_time=False
                )
            createOrderParameterPlot(
                in_file=args.in_file,
                fig_file=args.energy_fig_file.replace("energy","Q"),
                colour=cm.viridis,
                class_dict=class_dict,
                show=False,
                ax_lim='full',
                include_time=False,
                dr=0.1
                )
            createDensityPlot(
                in_file=args.in_file,
                fig_file=args.energy_fig_file.replace("energy","density"),
                colour=cm.viridis,
                class_dict=class_dict,
                show=False,
                ax_lim='full',
                include_time=False
                )
            for energy_type in ['Chaining','Hertzian','both']:
                try:
                    createEnergyPlot(
                        in_file=args.in_file,
                        fig_file=args.energy_fig_file,
                        energy_type=energy_type,
                        colour=cm.viridis,
                        class_dict=class_dict,
                        show=False,
                        bend_only=True,
                        ax_lim='full',
                        include_time=False
                        )
                except Exception as e:
                    print(e)

        except Exception as e:
            print(e)
            end=len(glob.glob(f"{args.base_path}/{args.data_dir}/biofilm_*.dat"))-1
            in_file=f"{args.base_path}/{args.data_dir}/biofilm_{end:05d}.dat"
            createFull2DPlot(
                in_file=in_file,
                fig_file=args.fig_file,
                colour=cm.viridis,
                class_dict=class_dict,
                show=False,
                ax_lim='full',
                include_time=False
            )
        quit()
    else:
        makeMovie(start=args.start_index,
                  end=args.end_index,
                  step=args.step,
                  base_path=args.base_path,
                  stem=args.stem,
                  data_dir=args.data_dir,
                  animation_dir=args.animation_dir,
                  max_cores=args.max_cores,
                  class_dict=class_dict,
                  track=True)
