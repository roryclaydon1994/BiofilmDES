# Standard libraries
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as linalg
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.patches import Arc

# Custom libraries

class SphericalBacterium(object):
    """docstring for SphericalBacterium."""

    sus_vis_radius_factor=1

    def __init__(self,cell_id,radius,
                      pos_x,pos_y,pos_z):
        self.cell_id          = int(cell_id)    # cell id
        self.radius           = radius          # spherocylinder radius
        self.diameter         = 2*radius        # spherocylinder diameter
        self.rcm = np.array([pos_x,pos_y,pos_z])
        self.colour = (5/255,106/255,117/255,0.8)

    def addElementToPlot(self,ax,colour=None,ax_rng=20,show_id=False):
        """
            Draw a spherical bacterium on an axis
            Parameters:
                ax: axis handle
                    The axis to draw the particle on
                colour: colour
                    Specify fill colour for this element
            Returns:
                perp_vector: float array, 3 elements
                    The vector perpendicular to the orientaion
        """

        if colour==None:
            colour=self.colour

        plot_radius=self.radius*SphericalBacterium.sus_vis_radius_factor

        # Represented as just a circle
        circle = plt.Circle(self.rcm,
                            radius=plot_radius,
                            fc=colour,edgecolor='k',
                            linewidth=5/ax_rng)
        ax.add_patch(circle)

        if show_id==True:
            ax.text(self.rcm[0],self.rcm[1],f"{self.cell_id}")

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        """
        :)
        """
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
