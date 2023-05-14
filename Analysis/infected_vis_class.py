"""
    Script to visualise infected bacteria from BiofilmDES

    This script defines the class Infected, child of RodShapedBacterium
    to represent the elements in the plot.

    Desctiption to follow
"""

# Add parent to path
import sys
sys.path.insert(0,'..')

# Custom modules
from RodShapedBacteria import RodShapedBacterium

class Infected(RodShapedBacterium):
    def __init__(self,cell_id,length,diameter,
                      com_vec_x,com_vec_y,com_vec_z,
                      orientation_x,orientation_y,orientation_z, lysis_period, burst_size,
                      virtual_center_x=0.0,virtual_center_y=0.0,virtual_center_z=0.0,
                      force_x=0.0,force_y=0.0,force_z=0.0,
                      torque_x=0.0,torque_y=0.0,torque_z=0.0,
                      neighbours=None,
                      projection="xy"):
        self.lysis_period = lysis_period
        self.burst_size = burst_size
        super(Infected, self).__init__(cell_id,length,diameter,
                          com_vec_x,com_vec_y,com_vec_z,
                          orientation_x,orientation_y,orientation_z,
                          virtual_center_x=0.0,virtual_center_y=0.0,virtual_center_z=0.0,
                          force_x=0.0,force_y=0.0,force_z=0.0,
                          torque_x=0.0,torque_y=0.0,torque_z=0.0,
                          neighbours=None,
                          projection="xy")
