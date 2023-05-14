"""
    Script to visualise phage from BiofilmDES

    This script defines the class Phage to represent the elements in the plot

    Desctiption to follow
"""

# Add parent to path
import sys
sys.path.insert(0,'..')

# Standard libraries
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as linalg
import pandas as pd
from FindDefects import intialiseElementsFromData
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

class Phage(object):
    """
        Control visualisation of phage

        Author: Rory Claydon

    """

    def __init__(self, cell_id,radius,pos_x,pos_y,pos_z):
        self.cell_id = cell_id
        self.rcm = np.array([pos_x,pos_y,pos_z])
        self.radius = 0.5
        self.colour = (181/255,4/255,66/255,0.8)

    def addElementToPlot(self,ax,colour=None,ax_rng=20.0,show_id=False):
        """
            Defines the representation of this object in 2d

            Parameters:
                ax: matplotlib axis
                    axis to add this phage to (expects equal aspect)
                colour: matplotlib colour
                    any type compatible with matplotlib colour
                ax_rng: float
                    max(x)-min(x) where x is the domain of ax
            Returns:
                mutated ax with this object added
        """
        if colour==None:
            colour=self.colour

        if show_id==True:
            ax.text(self.rcm[0],self.rcm[1],f"{self.cell_id}")
        circle = plt.Circle(self.rcm[:2],
                            radius=self.radius,
                            fc=colour,edgecolor='k',linewidth=5/ax_rng)
        ax.add_patch(circle)

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        patches=[]

        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        radius=0.5*height
        diameter=height
        tangent_axes=[0,1]
        colour=self.colour
        lw=0.1

        center=x0+0.5*width,y0+0.5*height

        circle = mpatches.Circle(center,
                                 radius=radius,
                                 fc=colour,edgecolor='k',lw=lw,
                                 transform=handlebox.get_transform())
        patches.append(circle)

        patch_collection = PatchCollection(patches,match_original=True)
        handlebox.add_artist(patch_collection)
        return patch_collection

def calcMSD(final_index=100):
    """
        Find the mean squared displacement of the phage

        Parameters:

        Returns:

    """
    phg_initial = intialiseElementsFromData(
        pd.read_csv("test_phage_motion_00000.dat",sep='\t'),
        Phage
        )

    phg_final = intialiseElementsFromData(
        pd.read_csv(f"test_phage_motion_{final_index:05d}.dat",sep='\t'),
        Phage
        )

    seps = np.zeros((len(phg_final),3))
    for ii,(phg_ii,phg_ff) in enumerate(zip(phg_initial,phg_final)):
        seps[ii] = phg_ff.rcm - phg_ii.rcm

    # Diffusion coefficient as defined in calcallforces for brownian motion
    n=2 #2D
    dt=0.1 #timestep used
    A=1
    D=A**2/2
    msd = np.sum(np.average(seps**2, axis=0))
    print(f"MSD is: {msd}") #Off by factor of 2?
    print(f"Expect: msd = 2nDt. Thus at this time and {n}D: {2*n*D*dt*final_index}")


if __name__ == "__main__":

    calcMSD()

    fig,ax = plt.subplots(1,1,figsize=[10,10])
    phg_data=pd.read_csv("test_phage_motion_00100.dat",sep='\t')
    all_phage = intialiseElementsFromData(phg_data,Phage)
    for phage in all_phage:
        phage.addElementToPlot(ax,ax_rng=40)
    ax.set_xlim([-20,20])
    ax.set_ylim([-20,20])
    plt.show()
