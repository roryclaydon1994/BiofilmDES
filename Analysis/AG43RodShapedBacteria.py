"""
    Chaining RodShapedBacterium element class

    To do:
        Remove recursion from the functions for finding chain head and tail.

    Author: Rory Claydon
"""

# Standard libraries
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg
import copy
from RodShapedBacteria import RodShapedBacterium

#third party libraries
from shapely.geometry import LineString

class AG43RodShapedBacterium(RodShapedBacterium):
    """
        Control visualisation of chaining susceptible in a biofilm

    """
    # cell_hash = {}

    def __init__(self,cell_id,length,radius,
                      pos_x,pos_y,pos_z,
                      ori_x,ori_y,ori_z,
                      springs,
                      virtual_center_x=0.0,
                      virtual_center_y=0.0,
                      virtual_center_z=0.0,
                      force_x=0.0,
                      force_y=0.0,
                      force_z=0.0,
                      torque_x=0.0,
                      torque_y=0.0,
                      torque_z=0.0,
                      neighbours=None):
        super().__init__(
            cell_id           = cell_id,
            length            = length,
            radius            = radius,
            pos_x             = pos_x,
            pos_y             = pos_y,
            pos_z             = pos_z,
            ori_x             = ori_x,
            ori_y             = ori_y,
            ori_z             = ori_z,
            virtual_center_x  = virtual_center_x,
            virtual_center_y  = virtual_center_y,
            virtual_center_z  = virtual_center_z,
            force_x           = force_x,
            force_y           = force_y,
            force_z           = force_z,
            torque_x          = torque_x,
            torque_y          = torque_y,
            torque_z          = torque_z,
            neighbours        = neighbours
        )
        self.springs = {}
        if springs != '[]':
            # print(springs)
            spring_split=springs.replace("[","").replace("]","").replace("(","").split('),')
            for ss in spring_split:
                if ss != '':
                    data = ss.split(',')
                    self.springs[int(data[0])] = [float(dd) for dd in data[1:]]

    def addElementToPlot(self,ax,colour=None,projection="xy",ax_rng=20,show_id=False):
        """
            Draw a susceptible on an axis, but now add spring links.
            Expects the hash list to already have been created.
            Parameters:
                ax: axis handle
                    The axis to draw the particle on
                colour: colour
                    Specify fill colour for this element
                projection: str
                    determines onto which plane to project the 3D biofilm
                    options: "xy","xz","yz"
                annotate: str
                    annotate element with this string
            Returns:
                perp_vector: float array, 3 elements
                    The vector perpendicular to the orientaion
        """
        super().addElementToPlot(ax,colour=colour,projection=projection,
                                 ax_rng=ax_rng,show_id=show_id)
        try:
            for id,params in self.springs.items():
                cell2 = self.cell_hash[id]
                p1=self.rcm  + params[0]*self.ori*self.length*0.5
                p2=cell2.rcm + params[1]*cell2.ori*cell2.length*0.5
                ov=np.array(params[2:])
                # This needs to be implemented once correct surface point selected
                # cv=p2-p1
                p1+=self.sus_vis_radius_factor*self.radius*ov
                p2-=self.sus_vis_radius_factor*cell2.radius*ov
                ax.plot([p1[0],p2[0]],[p1[1],p2[1]],'k',lw=5/ax_rng,zorder=3,
                        marker='o',mec='k',mfc='k',markersize=5/ax_rng
                        )
        except Exception as e:
            print(e)
            pass

    @staticmethod
    def createSpringHashList(cells):
        # AG43RodShapedBacterium.cell_hash={ cell.cell_id : cell for cell in cells }
        AG43RodShapedBacterium.makeHashList(cells)

    def getSpringGons(self,radiusfactor=1):
        if not AG43RodShapedBacterium.cell_hash:
            print("Error! Please create cell hash list first")
            quit()
        springs=[]
        for id,params in self.springs.items():
            cell2 = self.cell_hash[id]
            p1=self.rcm  + params[0]*self.ori*self.length*0.5
            p2=cell2.rcm + params[1]*cell2.ori*cell2.length*0.5
            ov=np.array(params[2:])
            # This needs to be implemented once correct surface point selected
            # cv=p2-p1
            p1+=radiusfactor*self.radius*ov
            p2-=radiusfactor*cell2.radius*ov
            springs.append(LineString([p1[:2],p2[:2]]))
        return springs
