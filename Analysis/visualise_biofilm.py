"""
    Create 2D projections of 3D outputs from BiofilmPhageDES

    This script contiains the plotting class and requires the files:
    RodShapedBacteria.py
    ChainingRodShapedBacteria.py
    which define the objects to be plotted.

    A quick example is provided at the bottom and will be run if this script is
    invoked directly. This will be removed soon.

    Expected usage:

    python3 visualise_biofilm.py sim_raw_data_dir file_stem save_dir

    Arguments:
        sim_raw_data_dir
            The path to the directory where the simulation data lives.

        file_stem
            I save files in the format filestem{%05d}.txt
            In other words, filestem is the leading part of the dumped sim data
            filename.

        save_dir
            directory to where you would like to direct anything saved in the
            script

    i.e. python3 visualise_biofilm.py data vis_biofilm_ output

    Intances of VisBiofilm need to be created from sim data in order to create
    the images.

    Functions (will be moved to staticmethod later) at the bottom create the
    images given instances of VisBiofilm.

    Plotting functions:

        createFull2DPlot
            creates a single projected 2D plot. Now has the capabilities to plot
            slices

        createLastFull2DPlot
            creates the final outputted 2D plot

        createBiofilmGrowingMovie
            generates all the plots required to make a movie.
            Currently uses all available cpus to parallelise over.

        createMovieFromFileIndices
            generate an image given the output timestep index

Author: Rory Claydon
"""

# Third party
# Experimental investigation of third party backend
# import mplcairo.base
# import mplcairo.qt

# Standard modules
import pandas as pd
import matplotlib as mpl
import re
import argparse
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.pyplot import axis, cm
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import glob
import os.path
import sys
import multiprocessing as mp
import pandas

# User defined classes
from RodShapedBacteria import RodShapedBacterium
from ChainingRodShapedBacteria import ChainingRodShapedBacterium
from FindDefects import Node

# User defined utility functions
from RodShapedBacteria import setProjection       # choose plane to visualise
from FindDefects import singleFrameFromCellData   # find defects in plane
from find_channels import isPointInCell           # check if point is in cell
from find_channels import getCavities             # identify cavities
from FindDefects import intialiseElementsFromData # create cell list from df
from FindDefects import cellPlaneIntersect        # check cell intersects with plane
try:
    from tylerUtilities import loadTylerData          # integrate LC hybrid solver data from Tyler
except Exception as e:
    pass

def setClassStaticMembers(log_file):
    hps=pd.read_csv(log_file,sep="\t").to_dict('list')
    if "Rd" in hps:
        RodShapedBacterium.depletion_radius=hps['Rd'][0]
    if "RodModE" in hps:
        RodShapedBacterium.elastic_moduli=hps['RodModE'][0]
        # print(f"Set {RodShapedBacterium.elastic_moduli=}")
    if "RodAspectRatio" in hps:
        RodShapedBacterium.smoothCutOff=hps['RodAspectRatio'][0]+1
        # print(f"Set {RodShapedBacterium.smoothCutOff=}")
    if "Kappa" in hps:
        ChainingRodShapedBacterium.kappa=hps['Kappa'][0]
        # print(f"Set {ChainingRodShapedBacterium.kappa=}")
    if "BendRig" in hps:
        ChainingRodShapedBacterium.bending_moduli=hps['BendRig'][0]
        # print(f"Set {ChainingRodShapedBacterium.bending_moduli=}")
    return hps

def colourByAngle(element,colour_scheme=cm.Greens):
    """
        colour elements by angle

        Parameters:
            element: instance of RodShapedBacterium or ChainingRodShapedBacterium class
                the element to be coloured in
            colour_scheme: callable (float) -> RGBA
                colour map to apply
        Returns:
            element_colour: RBGA
                The fill coulour for this element
    """
    theta = element.theta / np.pi

    # nematic so mod pi direction
    colour_index = int(np.floor(self.repeat_length*theta))
    return colour_scheme(colour_index)

def colourByHeight(element,min_z,max_z,colour_scheme=cm.Greens):
    """
        colour elements by z coordinate

        Parameters:
            element: instance of RodShapedBacterium or ChainingRodShapedBacterium class
                the element to be coloured in
            min_z (max_z): float
                minimum (maximum) z coordinate of all the cells
            colour_scheme: callable (float) -> RGBA
                colour map to apply
        Returns:
            element_colour: RBGA
                The fill coulour for this element
    """
    max_cap = 15
    max_z = max(max_z,max_cap)

    min_z = max(min_z,0)

    capped_colour=min(max(element.com_vec_z,min_z),max_z)-min_z
    colour = 1-np.log(1+capped_colour)/np.log(1+max_z-min_z)

    return colour_scheme(colour)

class VisBiofilm(object):
    """
        Visualisation of the biofilm class.

        Given a cell Class and an output file, visualise a projection onto the
        xy, xz or yz plane. Each cell Class has an add element method which
        defines how it will be visualised. Movies can also be created.

        Methods:

        Attributes:


    """

    def __init__(self,filename,Class,defect_filename="",repeat_length=8):
        """
            Visualise biofilms produced from data located in filename.

            Parameters:
                filename: str
                    name of the input file. Expects "\t" delimeted data.
                Class: object Class
                    Class of the data type to be plotted.
                defect_filename: str
                    name of the file containing the defect data
                repeat_length: int
                    Length of colour sequence (default is 8 for pastel scheme).
                    This is only a visual aid to coulour chains of bacteria.

        """
        self.repeat_length = repeat_length
        self.data = pd.read_csv(filename,sep="\t")
        self.defect_filename = defect_filename

        try:
            self.defect_data, self.defect_origin = loadTylerData(defect_filename)
        except Exception as e:
            print(e)
            print(f"failed loading defect: {defect_filename}")

        self.element_list = intialiseElementsFromData(self.data,Class)
        self.params={'font.family'           :'serif',
                     'text.usetex'           : True,
                     'axes.titlesize'        : 10,
                     'axes.labelsize'        : 10,
                     'xtick.labelsize'       : 10,
                     'ytick.labelsize'       : 10,
                     'legend.frameon'        : False,
                     'lines.linewidth'       : 1,
                     # 'legend.title_fontsize' : None,
                     'legend.fontsize'       : 10,
                     'font.size'             : 10,
                     "axes.formatter.useoffset":False,
                     "text.latex.preamble"   : [
                                                r"\usepackage[detect-all,locale=DE]{siunitx}",
                                                r"\usepackage[T1]{fontenc}",
                                                r"\usepackage{amsmath}",
                                               ]
                    }

    def setUpFigAx(self,fig_size,bk='k',projection=None):
        """
            Initialise fig and axes
            Parameters:
                fig_size: float
                    figure size in inces
                    see setFigSize function for conversion from pt to inces
        """
        self.fig,self.ax = plt.subplots(1,1,figsize=fig_size,
                                        constrained_layout=True,
                                        facecolor=bk,
                                        subplot_kw={'projection': projection})

    def setElementColourScheme(self,colour_scheme=cm.Pastel2):
        """
            This function needs updating

            Set the display colour of the elements
            Parameters:
                colour_scheme: colour map object from pyplot
                    The colour map used for display
        """
        evenly_spaced_interval = np.linspace(0, 1, self.repeat_length)
        colours = [colour_scheme(x) for x in evenly_spaced_interval]
        self.elements_colour_scheme = colour_scheme
        self.element_colours = colours
        return colours

    def setSaveName(self,save_name):
        """
            Set save name. Also checks if the path exists and if not creates it.
            Parameters:
                save_name: str
                    save name of the image created
        """
        filedir = os.path.dirname(save_name)
        # print(filedir)
        if filedir == '':
            pass
        elif not os.path.exists(filedir):
            os.makedirs(filedir)
        self.save_name = save_name

    def addElementsToPlot(self,projection="xy",ax_rng=20,colour_fn=colourByAngle):
        """
            Populate plot with cells.

            This is the default method for this purpose.
            The colour index used to be based on colouring cells according to
            which 3 degrees segment of the unit sphere they where aligned with.
            Currently there is no real meaning.

            Parameters:

                colour_fn: callable (cell_instance,colour_map) -> RGBA colour
                    determine how to colour elements
        """
        for ii,element in enumerate(self.element_list):
            element.addElementToPlot(self.ax,
                                     colour=colour_fn(element,self.elements_colour_scheme),
                                     projection=projection,
                                     ax_rng=ax_rng)

    def addChainedElementsToPlot(self,projection="xy",ax_rng=20.0,
                                 draw_springs=False,colour_chain=True):
        """
            Populate plot with chaining cells.

            Parameters:
                projection: str
                    determines onto which plane to project the 3D biofilm
                    options: "xy","xz","yz"
                ax_rng: float
                    the maximum range of any dimension of the plot. This is
                    required to scale the line width of the cell boundaries as
                    the domain increases
                draw_springs: bool
                    if True the spings are draw between all the cells which are
                    linked
                colour_chain: bool
                    if True then chains are all coloured the same. Otherwise
                    colour by height.

        """

        if colour_chain is True:
            # Adds cells to plot and colours them according to which chain they are in
            ChainingRodShapedBacterium.colourChains(self.ax,
                                             self.element_list,
                                             self.elements_colour_scheme,
                                             projection,
                                             ax_rng)
        else:
            max_z = max([cell.com_vec_z for cell in self.element_list])
            min_z = min([cell.com_vec_z for cell in self.element_list])
            colour_fn = lambda element,cm: colourByHeight(element,min_z,max_z,cm)
            self.addElementsToPlot(projection,ax_rng=ax_rng,colour_fn=colour_fn)

        # This draws the springs between the cells.
        if draw_springs is True:
            ChainingRodShapedBacterium.addLinksToElementsInPlot(
                self.ax,
                self.element_list
            )

    def addDefectFromFile(self,plane_loc=None,
                          eps=1e-4,projection="xy",ax_rng=20):
        """
            Add defects to plot

            Parameters:
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

            Returns:
                legend_elements: list of patches to add to legend

            Effect:
                Adds defects within the slice to the plot and colourises them based on charge
        """

        #switch selection depending on the projection being used
        # Roru: This will be changed to: setProjection to replace
        # select = []
        # axisSelect = 0
        # if projection == "xy":
        #     select = ['posX','posY']
        #     axisSelect = 'posZ'
        # elif projection == "xz":
        #     select = ['posX','posZ']
        #     axisSelect = 'posY'
        # elif projection == "yz":
        #     select = ['posY','posZ']
        #     axisSelect = 'posX'
        #
        # # Remove all defects outside of the appropriate plane
        # reducedDefectData = self.defect_data
        # if plane_loc is not None:
        #     reducedDefectData = self.defect_data[
        #         (abs(self.defect_data[axisSelect])-plane_loc > -eps)
        #         & (abs(self.defect_data[axisSelect])-plane_loc <  eps)
        #         ]
        #
        # # Get the appropriate axes from the selection
        # x = reducedDefectData[select[0]]
        # y = reducedDefectData[select[1]]
        # charge = reducedDefectData['charge']
        #
        # # Rory: size can be scaled by the axis size so the points always
        # # appear the same
        # print("trying to plot defects")
        # print(x)
        # self.ax.scatter(x,y,c=charge,s=20)
        # for x_curr,y_curr,charge_curr in zip(x,y,charge):
        #     if charge_curr >=0:
        #         colour="r"
        #     else:
        #         colour="b"
        #     circle = plt.Circle((x_curr,y_curr),
        #                         radius=5e-3*0.5*ax_rng,
        #                         fc=colour,
        #                         edgecolor='k',linewidth=5e-3*0.5)
        #     self.ax.add_patch(circle)
        # load Tyler's defect data
        charge, origin=loadTylerData(self.defect_filename)
        tyler_z_index=int(np.around(plane_loc-origin[0]))
        print(tyler_z_index)
        print(f"adding defects from {origin[0]+tyler_z_index} to plane at {plane_loc}")
        charge=charge[2]
        positive=np.where(charge>0)
        negative=np.where(charge<0)
        # vis_bio.ax.scatter(origin[2]+positive[2],origin[1]+positive[1],c='r')
        # vis_bio.ax.scatter(origin[2]+negative[2],origin[1]+negative[1],c='b')
        print(f"Found {len(positive[0])} +0.5 defects and {len(negative[0])} -0.5 defects")
        # print(f"The physical size of the gridspacing is: {ax_lim/Node.gridspacing}")
        #
        # # Try to reproduce on the fly
        # Node.epsilon_layer=eps

        # defect_container=singleFrameFromCellData(vis_bio.element_list,
        #                                          plane_loc=plane_loc)
        # print(f"Found {len(defect_container)} defects")
        #
        scal=15
        cmap = cm.get_cmap('bwr')
        ax_lim = self.getAxisLimits(self.element_list)
        legend_elements=[]
        for charge,d_pos in list(zip([-0.5,0.5],[negative,positive])):
            print(f"doing {charge}")
            d_colour=cmap(0.5*(1+charge))
            legend_elements.append(
                Line2D([0], [0], marker='o', color='w',
                       label=f'{charge:2.1f} defect',
                       markerfacecolor=d_colour,markeredgecolor='k', markersize=15)
                       )
            d_pos_arr=(origin[1:,np.newaxis]+np.asarray(d_pos)).T
            for y,x in d_pos_arr:
                # if z!=tyler_z_index:
                #     continue
                # defect_data = defect.getDefectData()
                # d_pos = defect_data[:3]
                # charge = defect_data[3]

                # Make a circle at the defect centre

                circle = plt.Circle((x,y),
                                    radius=0.1*ax_lim/scal,
                                    fc=d_colour,
                                    edgecolor='k',linewidth=1)
                self.ax.add_patch(circle)

                # Add a dashed line to indicate the size of the grid spacing
                # circle = plt.Circle((x,y),
                #                     radius=0.5*ax_lim/scal,
                #                     fc="none",
                #                     edgecolor=d_colour,linestyle="--",linewidth=2)
                # self.ax.add_patch(circle)
        print("Added defects")
        return legend_elements

    def findDefectFromElements(self,cell_list):
        """
            Find defect locations using a specific cell list

            Parameters:
            cell_list: list of Cell instances
                The cells for which to identify defects. It is assumed that these
                cells have already been correctly filtered into the appropriate
                layer with a preset tolerance.

        Effect:
            Adds defects within the slice to the plot and colourises them based on charge
        """

    @staticmethod
    def getAxisLimits(element_list,ax_init=20.0,ax_incr=2):
        """
            Resize axes to display all cells.

            Parameters:
                element_list: list of cells
                    cells to find the limits of
                ax_init: float
                    initial size of the axes
                ax_incr: float
                    Factor by which to increase the size of the axes
            Effect:
                sets the axis limits
        """

        max_position = max( [ max(element.rcm) for element in element_list ] )
        min_position = min( [ min(element.rcm) for element in element_list ] )
        print(f"min: {min_position} max: {max_position}")
        # use_position = max(max_position,ax_init)
        #
        # # Select the axis size as L = ax_init*ax_incr**inc_exp > max_pos
        # # with inc_exp the first large enough integer
        # inc_exp = int(
        #     np.ceil( np.log(use_position/ax_init) / np.log(ax_incr) )
        #     )
        #
        # ax_lim = ax_init * ( ax_incr ** inc_exp )

        ax_lim = 2*max(abs(max_position),abs(min_position)) + 5

        return ax_lim

def createPairInteractionFigure(vis_bio,format='pdf'):
        """
            Create a figure of the cells interacting.

            Parameters:
                visbio: VisBiofilm class instance
                    object which contains all the cell data
                format: str
                    file extension with which the image will be saved

        """
        assert len(vis_bio.element_list) == 2
        vis_bio.setSaveName(f"output/pair_interaction.{format}")

        vis_bio.setUpFigAx(fig_size = setFigSize(2*247))
        vis_bio.setElementColourScheme()
        for ii,element in enumerate(vis_bio.element_list):
            print(element)
            element.addElementToPlot(vis_bio.ax,colour=vis_bio.element_colours[ii])
            element.addHertzianInteractionToPlot(vis_bio.ax)

        vis_bio.ax.axis("scaled")
        vis_bio.ax.axes.get_xaxis().set_visible(False)
        vis_bio.ax.axes.get_yaxis().set_visible(False)
        vis_bio.ax.axis('off')
        vis_bio.fig.savefig(vis_bio.save_name,format=format,
                    bbox_inches='tight',
                    transparent=True)
        plt.show()

# This will be merged with the createFull2DPlot later, but should work for now
def createFull2DPlot(vis_bio,filename,plane_loc=None,eps=1e-4,projection="xy",
                     format='pdf',colour=cm.Greens,locate_channel=False):
    """
        Creates a plot of a slice through the biofilm centered on the chosen
        location and tolerance

        Parameters:
            vis_bio: VisBiofilm instance
                Contains the cell data and controls the appearence of the cells
            filename: str
                This is the filename which the image is saved to. (why not savename?)
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
            locate_channel: bool (default False)
                If true try to identify channels

        Effect:
            Saves the slice
    """
    vis_bio.setSaveName(f"{filename}.{format}")
    vis_bio.setUpFigAx(fig_size = setFigSize(2*247))
    vis_bio.setElementColourScheme(colour)

    tangent_axes, normal_axes = setProjection(projection)

    # Select cells to display which have an interesction wihthin eps of the
    # plane location
    if plane_loc is not None:
        vis_bio.element_list = [
            element for element in vis_bio.element_list
                    if cellPlaneIntersect(element,plane_loc,eps,normal_axes)
            ]

    # Sort cells so image appears as it would along the projected line of sight
    vis_bio.element_list = sorted(vis_bio.element_list,
                                  key=lambda cell: cell.rcm[normal_axes])

    # Choose axis size so that the axes don't change too often in movies.
    ax_lim = vis_bio.getAxisLimits(vis_bio.element_list)

    # The ax_rng is required to make the edge of the cells the right thickness
    vis_bio.addChainedElementsToPlot(projection,ax_rng=ax_lim,colour_chain=False)

    # Idenfity channel locations and add density plot to figure
    #### To do: move channel generation to a separate function
    if locate_channel is True:

        # add density heatmap
        # RodShapedBacterium.updateProjection(projection)
        # rho_field, rho_grid = RodShapedBacterium.computeDensity(vis_bio.element_list,dr=0.1)
        # CS=vis_bio.ax.contourf(*rho_grid,rho_field)
        # divider = make_axes_locatable(vis_bio.ax)
        # cax = divider.append_axes("right", size="5%", pad=0.05)
        # vis_bio.fig.colorbar(CS,cax=cax)

        vis_bio.save_name=vis_bio.save_name.strip(f".{format}")
        box_width=0.5
        threshold=1
        vis_bio.save_name+=f"|bw_{box_width}|th_{threshold}"
        vis_bio.save_name+=f".{format}"
        l,n,x,y=getCavities(vis_bio.element_list,
                            ax_rng=ax_lim,
                            box_width=box_width,
                            threshold=threshold)
        print(f"n: {n}")

        # hacked custom colour map
        viridis = cm.get_cmap('viridis', n)
        newcolors = viridis(np.linspace(0, 1, n))
        blank = np.array([0, 0, 0, 0])
        newcolors[0:2, :] = blank
        newcmp = ListedColormap(newcolors)

        # create plot
        c = vis_bio.ax.pcolormesh(y,x,l,cmap=newcmp,vmin=0,vmax=n,shading='auto')
        # divider = make_axes_locatable(vis_bio.ax)
        # cax = divider.append_axes("right", size="5%", pad=0.05)
        # vis_bio.fig.colorbar(c,cax=cax)

        # skip if no defect data generated
        # Should be a separate function
    try:
        legend_elements = vis_bio.addDefectFromFile(
            plane_loc,eps,projection,ax_rng=ax_lim
            )
        vis_bio.ax.legend(handles=legend_elements)
    except Exception as e:
        print(e)

    vis_bio.ax.axis("scaled")
    vis_bio.ax.axes.get_xaxis().set_visible(False)
    vis_bio.ax.axes.get_yaxis().set_visible(False)
    vis_bio.ax.axis('off')

    # projection in xy is approximately circular
    vis_bio.ax.set_xlim([-0.5*ax_lim,0.5*ax_lim])
    vis_bio.ax.set_ylim([-0.5*ax_lim,0.5*ax_lim])

    if projection != "xy":
        # min_z = min([cell.com_vec_z for cell in vis_bio.element_list])
        # max_z = max([cell.com_vec_z for cell in vis_bio.element_list])
        # vis_bio.ax.set_ylim([min_z,20*(1+max_z%20)])

        # Will think of a better way to do this later
        vis_bio.ax.set_ylim([None,None])

    vis_bio.fig.savefig(vis_bio.save_name,format=format,
                        bbox_inches='tight',
                        transparent=True,
                        dpi=300)
    plt.show()
    plt.close()

def setFigSize(fig_width_pt,ratio='Golden'):
    """
        Convert from pt to inches

        Parameters:
            fig_width_pt: float
                figure width in pt
            ratio:  str or float
                this is the aspect ratio. Use "Golden" for golden ratio or a
                float which gives: fig_height =fig_width*scale
        Returns:
            fig_size: list
                [fig_width,fig_height]
    """
    inches_per_pt = 1.0/72.27               # Convert pt to inches

    if ratio == 'Golden':
        scale = (np.sqrt(5)-1.0)/2.0        # Aesthetic ratio
    else:
        scale = ratio

    fig_width  =fig_width_pt*inches_per_pt  # width in inches
    fig_height =fig_width*scale             # height in inches
    # fig_height = fig_width*0.85           # height in inches
    fig_size = [fig_width,fig_height]
    # print("the fig size is: {}".format(fig_size))
    return fig_size

def createBiofilmGrowingMovie(sim_raw_data_dir = "data",
                              file_stem        = "vis_biofilm_",
                              save_dir         = 'output',
                              plot_class       = RodShapedBacterium,
                              addToExisting    = True,
                              projection       = "xy"):
    """
        Generate all the figures necessary to create a movie.

        This will try to use all cores available. Edit the NUM_PROCS variable
        to change this behaviour.

        Parameters:
            sim_raw_data_dir:
                The path to the directory where the simulation data lives.
            file_stem:
                I save files in the format filestem{%05d}.txt
                In other words, filestem is the leading part of the dumped sim data
                filename.
            save_dir:
                directory to where you would like to direct anything saved in the
                script
            plot_class: Class
                The type of element to be added to the plot:
                RodShapedBacterium or ChainingRodShapedBacterium
            addToExisting: bool
                If True, check if there is raw data in the sim_raw_data_dir which
                still needs to be visualised in output. Otherwise just start from
                0 again and overwrite the contents of output.
            projection: str
                determines onto which plane to project the 3D biofilm
                options: "xy","xz","yz"

        Effect:
            Makes movie figures
    """

    if addToExisting is True:
        # Find the number of images already created in the output directory and
        # take off 1 for 0 indexing
        pattern = f"{save_dir}/{file_stem}*.png"
        start_index=number_of_files = max(len(glob.glob(pattern)) - 1,0)
    else:
        start_index=0

    # Set up parallel image creation
    number_of_files = len(glob.glob(f"{sim_raw_data_dir}/vis_biofilm_*.txt"))

    # Create movie backwards as often these are the more intesting figures
    indices = [fileindex for fileindex in reversed(range(start_index,number_of_files))]

    # Set up threads and load balance
    print("Number of processors: ", mp.cpu_count())
    NUM_PROCS = mp.cpu_count()
    chunks = [indices[ii::NUM_PROCS] for ii in range(NUM_PROCS)]

    pool = mp.Pool(mp.cpu_count())

    results = pool.starmap(
        createMovieFromFileIndices,
        [(chunk,
          sim_raw_data_dir,
          file_stem,
          save_dir,
          plot_class,
          projection) for chunk in chunks]
        )

    pool.close()

def createBiofilmSliceMovie(filename,
                            defect_filename,
                            save_file_stem,
                            plot_class       = RodShapedBacterium,
                            projection       = "xy",
                            eps              = 1e-4,
                            format           = "png",
                            locate_channel   = False):
    """
        Generate all the figures necessary to create a movie of slicing through
        the biofilm in slices of set width

        This will try to use all cores available. Edit the NUM_PROCS variable
        to change this behaviour.

        Parameters:
            filename: str
                data file name
            save_file_stem: str
                save file name without file extension
            defect_filename: str
                name of the file containing the defect data
            plot_class: Class
                The type of element to be added to the plot:
                RodShapedBacterium or ChainingRodShapedBacterium
            projection: str
                determines onto which plane to project the 3D biofilm
                options: "xy","xz","yz"
            eps: float
                bacteria with centers of mass in |plane_loc - com[select]|<eps
                will be displayed where select picks the direction of projection
                e.g. if projection="xy" then select = 2 for the z direction
            locate_channel: bool (default False)
                If true try to identify channels

        Effect:
            Creates all the figures necessary to look through slices of the
            biofilm
    """
    vsb = VisBiofilm(filename,ChainingRodShapedBacterium,defect_filename)

    tangent_axes, normal_axes = setProjection(projection)

    if projection == "xy":
        col_name = "com_vec_z"
    elif projection == "xz":
        col_name = "com_vec_y"
    elif projection == "yz":
        col_name = "com_vec_x"
    else:
        print(f"projection {projection} not recognised")
        quit()

    # Get the correct com range to plot over
    h_max = vsb.data[col_name].max()
    h_min = vsb.data[col_name].min()
    print(f"range: [{h_min},{h_max}]")

    # Need to store this as the element list is mutated in plotting fn
    element_list = intialiseElementsFromData(vsb.data,plot_class)

    locs = np.arange(h_min,h_max,2*eps)
    for plane_loc in locs:
        print(f"plotting {plane_loc/h_max}")
        vsb.element_list = element_list
        savename = "{}|prj_{}|loc_{:3.3f}|eps_{:3.3f}".format(save_file_stem,
                                                              projection,
                                                              plane_loc,
                                                              eps)
        print(savename)
        createFull2DPlot(vsb,
                         projection=projection,
                         filename=savename,
                         format=format,
                         plane_loc=plane_loc,
                         eps=eps,
                         locate_channel=locate_channel)
        quit()

def createMovieFromFileIndices(indices,sim_raw_data_dir,file_stem,save_dir,
                                plot_class,projection,plane_loc=None,eps=1e-4):
    """
        Generate a figure using the time step index or indices

        Parameters:
            indices: int, array like
                the time step indices to create movies from
            sim_raw_data_dir:
                The path to the directory where the simulation data lives.
            file_stem:
                I save files in the format filestem{%05d}.txt
                In other words, filestem is the leading part of the dumped sim data
                filename.
            save_dir:
                directory to where you would like to direct anything saved in the
                script
            plot_class: Class
                The type of element to be added to the plot:
                RodShapedBacterium or ChainingRodShapedBacterium
            projection: str
                determines onto which plane to project the 3D biofilm
                options: "xy","xz","yz"
            plane_loc: float
                value in the projected direction to center the plots on
                e.g. if projection="xy" then this is the z direction
            eps: float
                bacteria with centers of mass in |plane_loc - com[select]|<eps
                will be displayed where select picks the direction of projection
                e.g. if projection="xy" then select = 2 for the z direction

        Effect:
            Makes the figures from the input files:
                sim_raw_data_dir/file_stem{indices}.txt
    """

    for fileindex in indices:
        # print( "progress: {:2.3f} \r".format( fileindex / number_of_files ),end="")

        filename = "{}/{}{:05d}.txt".format(sim_raw_data_dir,
                                            file_stem,
                                            fileindex)

        print(filename)

        vis_sim_biofilm = VisBiofilm(filename,plot_class)

        savename = "{}/{}{:05d}".format(save_dir,file_stem,fileindex)
        print(savename)

        createFull2DPlot(vis_sim_biofilm,
                         projection=projection,
                         filename=savename,
                         format='png',
                         plane_loc=plane_loc,
                         eps=eps)

def createLastFull2DPlot(sim_raw_data_dir = "data",
                         file_stem        = "vis_biofilm_",
                         save_dir         = 'output',
                         plot_class       = RodShapedBacterium,
                         projection       = "xy",
                         plane_loc        = None,
                         eps              = 1e-4):
    """
        Generate a figure from the last time step index

        Parameters:
            sim_raw_data_dir:
                The path to the directory where the simulation data lives.
            file_stem:
                I save files in the format filestem{%05d}.txt
                In other words, filestem is the leading part of the dumped sim data
                filename.
            save_dir:
                directory to where you would like to direct anything saved in the
                script
            plot_class: Class
                The type of element to be added to the plot:
                RodShapedBacterium or ChainingRodShapedBacterium
            projection: str
                determines onto which plane to project the 3D biofilm
                options: "xy","xz","yz"
            plane_loc: float
                value in the projected direction to center the plots on
                e.g. if projection="xy" then this is the z direction
            eps: float
                bacteria with centers of mass in |plane_loc - com[select]|<eps
                will be displayed where select picks the direction of projection
                e.g. if projection="xy" then select = 2 for the z direction

        Effect:
            Makes the figure from the input file:
                sim_raw_data_dir/file_stem{final_index}.txt
    """

    files = glob.glob(f"{sim_raw_data_dir}/vis_biofilm_*.txt")
    # print(re.findall("\d+","data/vis_biofilm_00103.txt"))
    # quit()
    indices = [ int(re.findall("\d+",filename)[-1]) for filename in files ]
    indices.sort()
    # number_of_files = len(glob.glob(f"{sim_raw_data_dir}/vis_biofilm_*.txt"))

    createMovieFromFileIndices([indices[-1]],
                               sim_raw_data_dir,
                               file_stem,
                               save_dir,
                               plot_class,
                               projection,
                               plane_loc=plane_loc,
                               eps=eps)

if __name__ == "__main__":

    # update to argparse
    sim_raw_data_dir = sys.argv[1]
    file_stem = sys.argv[2]
    save_dir = sys.argv[3]

    # Uncomment for example
    # filename = "data/vis_biofilm_00560.txt"
    # defectname = "data/defectout_00560.txt"
    # savename = "output/vis_biofilm_00560"
    # vis_sim_biofilm = VisBiofilm(filename,ChainingRodShapedBacterium,defectname)
    # createFull2DPlot(vis_sim_biofilm,plane_loc=0,projection="xy",eps=0.1,filename=savename,format='pdf')
    # createFull2DPlot(vis_sim_biofilm,projection="xz",filename=savename,format='pdf')

    # createLastFull2DPlot(sim_raw_data_dir,
    #                      file_stem,
    #                      save_dir,
    #                      ChainingRodShapedBacterium,
    #                      projection="xy",
    #                      plane_loc=0.0,
    #                      eps=1e-5)

    # ---------------------- Set hyperparameters -------------------------------
    aLAYERHEIGHT          = 1
    aGRIDSPACING          = 20
    aSEARCHSPACING        = 0.1
    EPSILON_LAYER         = 1
    EPSILON_LOOP          = 0.5
    EPSILON_TOPOTHRESHOLD = 0.2
    R_NEIGHBOUR           = 8
    R_DENSITY             = 0.5*R_NEIGHBOUR
    DELTA_DENSITY         = 0
    N_LOOPMIN             = 1

    # Set Node with these I/O and defect calculation hyperparameter data
    Node.setHyperParams(aLAYERHEIGHT,
                        aGRIDSPACING,
                        aSEARCHSPACING,
                        epsilon_layer=EPSILON_LAYER,
                        epsilon_loop=EPSILON_LOOP,
                        epsilon_topothreshold=EPSILON_TOPOTHRESHOLD,
                        r_neighbour=R_NEIGHBOUR,
                        r_density=R_DENSITY,
                        delta_density=DELTA_DENSITY,
                        n_loopmin=N_LOOPMIN
                        )

    index=300
    createBiofilmSliceMovie("data/vis_biofilm_{:05d}.txt".format(index),
                            "tyler_output/charg_biofilm_{:05d}.txt".format(index),
                            "output/movie_{0:05d}/meeting_biofilm_{0:05d}".format(index),
                            ChainingRodShapedBacterium,
                            projection="xy",
                            eps=0.1,
                            locate_channel=True)

    # createBiofilmGrowingMovie(sim_raw_data_dir,file_stem,
    #                           save_dir,ChainingRodShapedBacterium,False)
