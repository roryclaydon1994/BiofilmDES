# Standard modules
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import pandas as pd
import os
import sys

# Third party
from shapely.geometry import Polygon

# Custom modules
from RodShapedBacteria import RodShapedBacterium
from findCriticalLengthProb import getCells
import utilities as ut
import fastCellPlotting as fcp
from DistributionFunctions import computeColonyContour, getColonyDensity
from generalPlotting import addNematicDirector

ut.setMPL()
save_dir = "/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/figures/overlaps_gpu/"


def plot_overlaps(par_overlap_file, non_par_overlap_file):
    RodShapedBacterium.sus_vis_radius_factor = 0.7

    for lbl, file in [['par', par_overlap_file], ['non', non_par_overlap_file]]:
        cells = getCells(file)

        fig, ax = plt.subplots(1, 1, figsize=ut.getFigsize(cols=2),
                               constrained_layout=True,
                               facecolor='k')

        dat = pd.read_csv(file, sep='\t')
        maxx = dat['pos_x'].max() + 3
        minx = dat['pos_x'].min() - 3
        maxy = dat['pos_y'].max() + 3
        miny = dat['pos_y'].min() - 3
        fcp.addAllCellsToPlot(cells, ax, ax_rng=maxx - minx, show_id=True, ec='w')

        ax.axis('scaled')
        ax.axis('off')
        ax.set_xlim([minx, maxx])
        ax.set_ylim([miny, maxy])
        stem = f"{save_dir}/plotCells_{file.replace('/', '_')}"
        stem.strip('.dat')

        fig.savefig(f"{stem}_{lbl=}.pdf",
                    transparent=False,
                    bbox_inches='tight'
                    )
        plt.show()

    fig, ax = plt.subplots(1, 1, figsize=ut.getFigsize(cols=2),
                           constrained_layout=True)
    both_overlaps = []
    for lbl, file in [['par', par_overlap_file], ['non', non_par_overlap_file]]:
        cells = getCells(file)
        overlaps = RodShapedBacterium.getOverlaps(cells)
        both_overlaps.append(overlaps)
        print(f"{lbl=} {np.max(overlaps)=}")
        # data, bins = np.histogram(overlaps, bins=np.linspace(0, 0.5, 20), density=True)
        # plt_dict[lbl] = {'bins': bins, 'data': data}
    ax.hist(both_overlaps, bins=10, density=True, label=['better_parallel', 'basic_parallel'], alpha=0.8)
    plt.legend(loc='upper right')
    stem = f"{save_dir}/plotCells_{file.replace('/', '_')}"
    stem.strip('.dat')

    fig.savefig(f"{stem}_overlaps.pdf",
                transparent=False,
                bbox_inches='tight'
                )
    plt.show()


if __name__ == "__main__":

    par = "/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/SimOutput/ChainingBucklingTransition/run17/repeat1" \
          "/final_00107.dat"
    non = "/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/SimOutput/ChainingBucklingTransition/run21/repeat1" \
          "/final_00103.dat"
    os.makedirs(save_dir, exist_ok=True)
    plot_overlaps(par, non)


    # ff = 'final_00200.dat'
    # par = f"/home/rory/PhD_sync/CUDA/cell_linked_list/par_output/run1/{ff}"
    # non = f"/home/rory/PhD_sync/CUDA/cell_linked_list/non_par_output/run1/{ff}"
    # os.makedirs(save_dir, exist_ok=True)
    # plot_overlaps(par, non)

    # ff = 'final_00200.dat'
    # par = f"/home/rory/PhD_sync/CUDA/cell_linked_list/GeneratedOutput/output_new_parallel/{ff}"
    # non = f"/home/rory/PhD_sync/CUDA/cell_linked_list/GeneratedOutput/output_old_parallel/{ff}"
    # os.makedirs(save_dir, exist_ok=True)
    # plot_overlaps(par, non)
