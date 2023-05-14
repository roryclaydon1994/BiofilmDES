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
from matplotlib.collections import PatchCollection
import numpy as np
import numpy.linalg as linalg
import pandas as pd
import argparse
import os
import subprocess
import multiprocessing as mp
from joblib import Parallel, delayed
from copy import copy
import re
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
fontprops = fm.FontProperties(size=12)

# Third party
from hidden_prints import HiddenPrints

# Custom modules
from FindDefects import intialiseElementsFromData
from visualise_infected_biofilm import VisInfectedBiofilm
from visualise_biofilm import setFigSize
from phage_vis_class import Phage
from infected_vis_class import Infected
from RodShapedBacteria import RodShapedBacterium
from RodShapedBacteria import setProjection  # choose plane to visualise

import matplotlib as mpl
params={'font.family'           :'serif',
        'text.usetex'           : True,
        'axes.titlesize'        : 10,
        'axes.labelsize'        : 10,
        'xtick.labelsize'       : 10,
        'ytick.labelsize'       : 10,
        'legend.frameon'        : False,
        'lines.linewidth'       : 1.5,
        'legend.handlelength'   : 1.5,
        'legend.columnspacing'  : 0.75,
        'legend.handletextpad'  : 0.5,
        # 'legend.title_fontsize' : None,
        'legend.fontsize'       : 10,
        'font.size'             : 10,
        'legend.columnspacing'  : 1.0,
        "axes.formatter.useoffset":False,
        "text.latex.preamble"   : [r"\usepackage[detect-all,locale=DE]{siunitx}",r"\usepackage[T1]{fontenc}",r"\usepackage{amsmath}"]
       }
mpl.rcParams.update(params)

class AnyObject(object):
    pass

def makeMovie(start, end, step=1, max_cores=mp.cpu_count(),base_path="",data_dir="data",animation_dir="animation",track=True):
    # Convert to parallel
    # for ii in range(args.input_index,args.upper_index+1):
    #     args.sus_file = f"data/main_bacteria_{ii:05d}.dat"
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
    sus_data = pd.read_csv(f"{base_path}/{data_dir}/main_bacteria_{0:05d}.dat", sep='\t')
    all_sus = intialiseElementsFromData(sus_data, RodShapedBacterium)
    ax_lim = np.max(VisInfectedBiofilm.getAxisLimits(all_sus)) + 5
    print(ax_lim)

    if track==True:
        ax_lim=None
    else:
        ax_lim=512

    """
        Generate list of tuples for each input to createFull2DPlot in the
        order
        sus_file,inf_file,phg_file,fig_file,sus_class=RodShapedBacterium,
        inf_class=Infected,phg_class=Phage,plane_loc=None,eps=1e-4,
        projection="xy",colour=cm.Greens,show=False
    """
    par_args_list = [
        {
            "sus_file": f"{base_path}/{data_dir}/main_bacteria_{ii:05d}.dat",
            "inf_file": f"{base_path}/{data_dir}/main_infected_{ii:05d}.dat",
            "phg_file": f"{base_path}/{data_dir}/main_phage_{ii:05d}.dat",
            "fig_file": f"{base_path}/{animation_dir}/infected_biofilm_{ii:05d}.png",
            "colour": cm.viridis,
            "show": False,
            "ax_lim": ax_lim,
            "dpi": 300
        } for ii in range(start, end + 1, step)]
    # Convert to Parallel
    print(f"Running with {max_cores} cores")
    # Parallel(n_jobs=max_cores, verbose=10, backend="multiprocessing")(
    #     delayed(
    #         suppressOutput(createFull2DPlot))(**par_args)
    #     for par_args in par_args_list
    # )
    Parallel(n_jobs=max_cores, verbose=10, backend="multiprocessing")(
        delayed(createFull2DPlot)(**par_args)
        for par_args in par_args_list
    )
    # with mp.Pool(processes=args.max_cores) as pool:
    #     results = pool.starmap(
    #         createFull2DPlot,
    #         [par_args for par_args in par_args_list]
    #         )

    if step == 1:
        movie_command = (
                f"ffmpeg -y -r 10 -i {base_path}/{animation_dir}/infected_biofilm_%05d.png"
                + " -c:v libx264 -vf fps=25 -pix_fmt yuv420p"
                + f" {base_path}/{animation_dir}/infected_biofilm.mp4"
        )
    else:
        movie_command = (
                "ffmpeg -y -f image2 -pattern_type glob -r 10"
                + f" -i '{base_path}/{animation_dir}/infected_biofilm_*.png'"
                + " -c:v libx264 -vf fps=25 -pix_fmt yuv420p"
                + f" {base_path}/{animation_dir}/infected_biofilm.mp4"
        )
    os.system(movie_command)


def suppressOutput(func):
    def noPrint(*args, **kwargs):
        with HiddenPrints():
            out = func(*args, **kwargs)
        return out

    return noPrint


def createFull2DPlot(sus_file, inf_file, phg_file, fig_file, sus_class=RodShapedBacterium,
                     inf_class=Infected, phg_class=Phage, plane_loc=None, eps=1e-4,
                     projection="xy", colour=cm.Greens, show=False, dpi=100,
                     ax_lim=None):
    """
        Creates a plot of a slice through the biofilm centered on the chosen
        location, within tolerance (eps) in the given plane (projection)

        Parameters:
            sus_file: str
                Load susceptible data
            inf_file: str
                Load infected data
            phg_file: str
                Load phage data
            sus_class, inf_class, phg_class: visualisation classes
                These define how objects will appear on the axes
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

        Effect:
            Saves the slice
    """

    # Set up axes and load susceptible
    load_dict = {
        "RodShapedBacterium": [sus_file, sus_class],
        "Infected": [inf_file, inf_class],
        "Phage": [phg_file, phg_class]
    }
    vis_bio = VisInfectedBiofilm(load_dict)

    vis_bio.setSaveName(fig_file)
    vis_bio.setUpFigAx(fig_size=setFigSize(2 * 247))

    # If left blank, use class default colour schemes
    vis_bio.setElementColourScheme()

    # ------------------------ Quasi 2D considerations -------------------------
    tangent_axes, normal_axes = setProjection(projection)

    # Select cells to display which have an interesction wihthin eps of the
    # plane location
    if plane_loc is not None:
        for element_list in vis_bio.species_dict.values():
            element_list = [
                element for element in element_list
                if cellPlaneIntersect(element, plane_loc, eps, normal_axes)
            ]

    # Sort cells so image appears as it would along the projected line of sight
    for element_list in vis_bio.species_dict.values():
        element_list = sorted(element_list,
                              key=lambda cell: cell.rcm[normal_axes])

    # ---------------------------------------------------------------------------

    # Choose axis size so that all elements are displayed
    if ax_lim == None:
        ax_lim = 0
        for element_list in vis_bio.species_dict.values():
            try:
                ax_lim = max(vis_bio.getAxisLimits(element_list), ax_lim)
            except Exception as e:
                print(e)
                for key, value in vis_bio.species_dict.items():
                    print(f"{key} has {len(value)} elements")
        ax_lim += 5
    ax_lim=2*(int(ax_lim/2)+1)
    ax_lim=min(ax_lim,200)
    # The ax_rng is required to make the edge of the cells the right thickness
    for ii, element_list in enumerate(vis_bio.species_dict.values()):
        colour_int = 150 * ii
        print(colour_int)
        for cell in element_list:
            # print(cell.ori)
            if ( cell.cell_id==334 ):
                colour_int = 450
            else:
                colour_int = 150 * ii
            cell.addElementToPlot(vis_bio.ax, ax_rng=ax_lim,
                                  colour=colour(colour_int))
            cell.colour=colour(colour_int)

    # Create a nice legend
    try:
        sus_plot_cell=copy(vis_bio.species_dict['RodShapedBacterium'][0])
        inf_plot_cell=copy(vis_bio.species_dict['Infected'][0])
        legend_labels= [f"RodShapedBacterium : {len(vis_bio.species_dict['RodShapedBacterium'])}",
                        f"Infected :    {len(vis_bio.species_dict['Infected'])}"]
    except Exception as e:
        sus_plot_cell=RodShapedBacterium(0,0,0,0,0,0,1,0,0)
        sus_plot_cell.colour=colour(0)
        inf_plot_cell=Infected(0,0,0,0,0,0,1,0,0,0,0)
        inf_plot_cell.colour=colour(150)
        legend_labels= [f"RodShapedBacterium : {len(vis_bio.species_dict['RodShapedBacterium'])}",
                        f"Infected :    {0}"]
    legend_elements=[sus_plot_cell,
                     inf_plot_cell]

    try:
        phg_plot_cell=copy(vis_bio.species_dict['Phage'][0])
        legend_elements.append(phg_plot_cell)
        legend_labels.append(f"Phage : {len(vis_bio.species_dict['Phage'])}")
    except Exception as e:
        phg_plot_cell=Phage(0,0,0)
        phg_plot_cell.colour=colour(300)
        legend_elements.append(phg_plot_cell)
        legend_labels.append(f"Phage : {0}")

    handler_map={RodShapedBacterium: sus_plot_cell,
                 Infected: inf_plot_cell,
                 Phage: phg_plot_cell}
    try:
        time=int(re.findall(r'\d+',os.path.split(sus_file)[1])[0])*0.05
        title=f"Time: {time:2.2f} hrs"
    except Exception as e:
        print(e)
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
    scalebar = AnchoredSizeBar(vis_bio.ax.transData,
                           10,
                           r'\SI{10}{\micro\meter}',
                           'lower right',
                           pad=0.1,
                           color='black',
                           frameon=False,
                           size_vertical=1,
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
                        transparent=True,
                        dpi=dpi)
    if show == True:
        plt.show()
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='visualise_phage.py',
        usage='%(prog)s [options] path',
        description='Visualisation script for infected biofilms',
        epilog='Authors: Ricardo Del Rio, Rory Claydon'
    )
    parser.add_argument("-I", "--input_index", type=int, required=True,
                        help="load data from"
                             + " data/main_{species}_{input_index:05d}.dat")
    parser.add_argument("-U", "--upper_index", type=int, required=False,
                        default=None,
                        help="load data from"
                             + " data/main_{species}_{input_index:05d}.dat")
    parser.add_argument("--step", type=int, required=False,
                        default=1,
                        help="skip step number outputs")
    parser.add_argument("-S", "--sus_file", type=str, required=False, default=None,
                        help="susceptible data filename")
    parser.add_argument("-Q", "--inf_file", type=str, required=False, default=None,
                        help="infected data filename")
    parser.add_argument("-P", "--phg_file", type=str, required=False, default=None,
                        help="phage data filename")
    parser.add_argument("-F", "--fig_file", type=str, default=None,
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
    parser.add_argument("-BD", "--base-dir", type=str, required=False, default="",
                        action="store", dest="base_path",
                        help="default directory to look for inputs")
    parser.add_argument("-DD", "--data-dir", type=str, required=False, default="data",
                        action="store", dest="data_dir",
                        help="default directory name in base for inputs")
    parser.add_argument("-AD", "--anim-dir", type=str, required=False, default="animation",
                        action="store", dest="animation_dir",
                        help="default directory name in base for inputs")
    parser.add_argument("--fig-format", type=str, required=False, default="png",
                        action="store", dest="fig_format",
                        help="set figure format")
    parser.add_argument("--sus-vis-radius-factor", type=float, required=False, default=1,
                        action="store", dest="sus_vis_radius_factor",
                        help="Set visualisation width this fraction of real width for bacteria")
    args = parser.parse_args()

    RodShapedBacterium.sus_vis_radius_factor=args.sus_vis_radius_factor

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

    if args.upper_index == None:
        if args.sus_file == None:
            args.sus_file = f"{args.base_path}/{args.data_dir}/main_bacteria_{args.input_index:05d}.dat"

        if args.inf_file == None:
            args.inf_file = f"{args.base_path}/{args.data_dir}/main_infected_{args.input_index:05d}.dat"

        if args.phg_file == None:
            args.phg_file = f"{args.base_path}/{args.data_dir}/main_phage_{args.input_index:05d}.dat"

        if args.fig_file == None:
            args.fig_file = f"{args.base_path}/{args.animation_dir}/infected_biofilm_{args.input_index:05d}.{args.fig_format}"

        createFull2DPlot(
            sus_file=args.sus_file, inf_file=args.inf_file,
            phg_file=args.phg_file, fig_file=args.fig_file,
            colour=cm.viridis, show=True
        )
    else:
        makeMovie(start=args.input_index,
                  end=args.upper_index,
                  step=args.step,
                  base_path=args.base_path,
                  data_dir=args.data_dir,
                  animation_dir=args.animation_dir,
                  max_cores=args.max_cores)
