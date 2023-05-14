"""
    test script for FindDefects.py and FindDefectsNeighbours.py

    functions:
        runTests
            tests a reference configuration with an m charged defect.
            the tested function is singleFrameFromCellData from FindDefects.py

    Direct invocation sets some basic hyperparameter choices used for these
    tests. The number of successes if reported. then assertions check which if
    any of the tests fails.

"""
# Standard libraries
import sys
import math
import numpy as np
import numpy.linalg
import pandas as pd
from typing import TypeVar, Type
import re
import numpy.random
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.optimize as scopt

# User defined
from RodShapedBacteria import RodShapedBacterium
from ChainingRodShapedBacteria import ChainingRodShapedBacterium
from visualise_biofilm import intialiseElementsFromData, cellPlaneIntersect
from FindDefects import Node,singleFrameFromCellData

def runTest(m=0.5,hyperparameter_dict=None,num_trials=5):
    """
        Check we can find a defect with charge m

        Parameters:
        m: float
            defect charge
        hyperparameter_dict: dictionary
            Set the hyperparameters in the Node class for defect finding
            calculation
        num_trials: int
            number of trials to try to find the best defect location

        Effect:
            Successful completion implies that the minus half defect was correctly
            identified. Assertions check the charge is correct too.
    """
    if hyperparameter_dict is not None:
        # User would like to change the hyperparameters
        Node.setHyperParams(**hyperparameter_dict)

    particle_list = []
    radii = [2,6]   # Set a radius to place the particles on

    ori_theta = lambda theta: np.array([np.cos(m*theta),np.sin(m*theta),0])

    for ii,radius in enumerate(radii):
        thetas = np.linspace(0,2*np.pi,(ii+1)*15,endpoint=False)
        for theta in thetas:
            length=np.random.uniform(2,3)
            ori_pert=np.random.uniform(-0.1*np.pi,0.1*np.pi,2)
            com_pert=np.random.uniform(-0.5,0.5,3)

            rcm=radius*np.array([np.cos(theta),np.sin(theta),0])
            com_vec_x,com_vec_y,com_vec_z=rcm+com_pert

            ori = ori_theta(theta)
            ori = ori / np.linalg.norm(ori)
            orientation_x,orientation_y,orientation_z=ori

            particle_list.append(RodShapedBacterium(cell_id=0,
                                             length=length,
                                             radius=0.5,
                                             pos_x=com_vec_x,
                                             pos_y=com_vec_y,
                                             pos_z=com_vec_z,
                                             ori_x=orientation_x,
                                             ori_y=orientation_y,
                                             ori_z=orientation_z))

    find_defect = lambda r: -abs(RodShapedBacterium.computeLocalCharge(particle_list,r,dr=2**(-7)))
    limits=np.array([[-np.max(radii),np.max(radii)],
                     [-np.max(radii),np.max(radii)]
                     ])
    best_rs = np.zeros((num_trials,2))
    for trial in range(num_trials):
        min_fun = 1e20
        for _ in range(20):
            # Create initial guess
            x0_x = np.random.uniform(*limits[0])
            x0_y = np.random.uniform(*limits[1])
            x0 = [x0_x,x0_y]

            # Find the minima from this guess
            res = scopt.minimize(find_defect,x0)

            if res.fun < min_fun:
                # Save the position if this was the best guess from this set
                min_fun = res.fun
                best_rs[trial] = res.x

    # best_r = r_maxima.mean(axis=0)
    # print(best_r)
    best_r = best_rs.mean(axis=0)
    std_r  = best_rs.std(axis=0)
    print(f"r={best_r}$\pm${std_r}")
    # Check the cell configuration and expected defect location
    fig,ax=plt.subplots(1,1,figsize=[10,10])
    ax.axis("scaled")
    ax_rng = 2*max(radii) + 4
    ax.set_xlim([-0.5*ax_rng,0.5*ax_rng])
    ax.set_ylim([-0.5*ax_rng,0.5*ax_rng])

    q_field, Q, q_grid = RodShapedBacterium.computeChargeDensity(particle_list,dr=0.1)
    CS=ax.contourf(*q_grid,q_field)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(CS,cax=cax)

    for particle in particle_list:
        # print(particle,particle.ori)
        # ax.scatter(*particle.rcm[0:2])
        # lower=particle.rcm[0:2]-0.5*particle.ori[0:2]*particle.length
        # upper=particle.rcm[0:2]+0.5*particle.ori[0:2]*particle.length
        # ax.plot([lower[0],upper[0]],[lower[1],upper[1]])
        particle.addElementToPlot(ax,ax_rng=ax_rng)

    # Add the expected defect location
    circle = plt.Circle((0,0),
                        radius=0.25,
                        fc='m',
                        edgecolor='k',linewidth=1)
    ax.add_patch(circle)

    # Add the predicted defect location
    cmap = cm.get_cmap('bwr')
    q_defect = RodShapedBacterium.computeLocalCharge(particle_list,best_r,dr=2**(-7))
    color = cmap(0.5*(1+np.sign(q_defect)))
    circle = plt.Circle(best_r,
                        radius=0.25,
                        fc=color,
                        edgecolor='k',linewidth=1)
    ax.add_patch(circle)

    circle = plt.Circle(best_r,
                        radius=np.linalg.norm(std_r),
                        fc="none",
                        edgecolor=color,linestyle="-.",linewidth=2)
    ax.add_patch(circle)

    ax.annotate(xy=(best_r[0],best_r[1]),
                xytext=(0.4,0.4),
                textcoords='figure fraction',
                s="q={:2.2f}".format(q_defect),
                color='w',
                fontsize=12,
                horizontalalignment='left',verticalalignment='top',
                bbox=dict(facecolor=color, alpha=0.5),
                arrowprops=dict(facecolor=color, shrink=0.05))

    ax.text(0.01,0.01,
            s="r=[{:2.1f},{:2.1f}]$\pm$[{:2.1f},{:2.1f}]".format(*best_r,*std_r),
            color='w',
            fontsize=12,
            horizontalalignment='left',verticalalignment='bottom',
            transform=ax.transAxes,
            bbox=dict(facecolor=color, alpha=0.5))

    # Find the defect from this cell configuration
    defect_container=singleFrameFromCellData(particle_list,plane_loc=0)
    for defect in defect_container:
        defect_data = defect.getDefectData()
        d_pos = defect_data[:3]
        charge = defect_data[3]
        d_colour=cmap(0.5*(1+charge))
        circle = plt.Circle((d_pos[0],d_pos[1]),
                            radius=0.25,
                            fc=d_colour,
                            edgecolor='k',linewidth=1)
        ax.add_patch(circle)
        circle = plt.Circle((d_pos[0],d_pos[1]),
                            radius=0.5*ax_rng/Node.gridspacing,
                            fc="none",
                            edgecolor=d_colour,linestyle="-.",linewidth=2)
        ax.add_patch(circle)
        ax.text(0.99,0.01,
                s="r=[{:2.1f},{:2.1f}]$\pm$[{:2.1f},{:2.1f}]".format(d_pos[0],
                                                                     d_pos[1],
                                                                     ax_rng/Node.gridspacing,
                                                                     ax_rng/Node.gridspacing),
                color='w',
                fontsize=12,
                horizontalalignment='right',verticalalignment='bottom',
                transform=ax.transAxes,
                bbox=dict(facecolor=d_colour,alpha=0.5))
        ax.annotate(xy=(d_pos[0],d_pos[1]),
                    xytext=(0.6,0.6),
                    textcoords='figure fraction',
                    s="q={:2.2f}".format(charge),
                    color='w',
                    fontsize=12,
                    horizontalalignment='right',verticalalignment='top',
                    bbox=dict(facecolor=d_colour, alpha=0.5),
                    arrowprops=dict(facecolor=d_colour, shrink=0.05))

    fig.savefig(f"defect_test_{m}.png",format="png",bbox_inches='tight')
    # plt.show()
    plt.close()

    return best_r,q_defect

if __name__ == "__main__":

    # ---------------------- Set hyperparameters -------------------------------
    aLAYERHEIGHT          = 1
    aGRIDSPACING          = 6
    aSEARCHSPACING        = 0.1
    EPSILON_LAYER         = 1
    EPSILON_LOOP          = 0.5
    EPSILON_TOPOTHRESHOLD = 0.2
    R_NEIGHBOUR           = 5
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
    # ------------------------- Run the tests ----------------------------------
    outputs=[]
    for m in [-0.5,0.5]:
        outputs.append(runTest(m=m))
        print(f"For charge {m} at [{0},{0}] we calculate {outputs[-1][1]} at {outputs[-1][0]}")
        assert np.sign(outputs[-1][1])==np.sign(m), f"Charge has wrong sign"

    # print(outputs)
    # successes = [abs(e[0]-e[1])<1e-3 for e in outputs]
    # for test in range(len(outputs)):
    #     assert successes[test]==True, f"Failed for charge {outputs[test]}"
    # print("Success!")
