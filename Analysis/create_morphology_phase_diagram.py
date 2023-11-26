"""

"""

# Standard modules
import os
from typing import Iterable, Optional
import functools
import glob

# Third party modules
import numpy as np
from nptyping import NDArray
import matplotlib.pyplot as plt
import pandas as pd

# Custom modules
import utilities as ut
from fastCellPlotting import addAllCellsToPlot
from DistributionFunctions import computeColonyContour
import RodShapedBacteria

def load_cells_from_file_generator(files: Iterable[str], delimeter: str = '\t'):
    """
    Load cell data as numpy array.

    At present, assumes the file is a plain txt file/csv with the first row the column titles:
    type, id, length, radius, pos_x, pos_y, pos_z, ori_x, ori_y, ori_z, lower_link, upper_link

    :param files: 
    :param delimeter: delimeter used in the csv file
    :returns: cell data as numpy array
    """
    for file in files:
        cell_data = np.loadtxt(
            file, skiprows=1, delimiter=delimeter, usecols=range(2, 10)
        )
        yield cell_data


class MorphologyPhaseDiagram:
    """
    Handle creation of a 2D morphology phase diagram.
    """

    def __init__(self, x, y, files, labels: Optional[Iterable]) -> None:
        """
        :param x: x-axis independent variable, assumed ascending order if not categorical
        :param y: y-axis independent variable, assumed ascending order if not categorical
        :param files:
            files corresponding to each x, y in row major order, i.e.
            if there are n_x, n_y independent variables for x,y respectively, then the phase space
            coordinates of files at the following indices are given by
            file index | x index | y index
            ------------------------------
                 0     |    0    |    0
               n_x-1   |  n_x-1  |    0
                n_x    |    0    |    1
            and so on.

        """
        self.x: Iterable = x
        self.y: Iterable = y
        self.files: Iterable[str] = files
        self.labels: Iterable[str] = labels

        self.n_x = len(self.x)
        self.n_y = len(self.y)
        if len(files) != self.n_x * self.n_y:
            raise ValueError("Not enourgh files for number of independent paramerters")

        self.box_widths = np.zeros(self.n_x)
        self.box_heights = np.zeros(self.n_y)

    def get_flat_index_from_row_col(self, x_index, y_index):
        return x_index + self.n_x * y_index

    def get_row_col_from_flat_index(self, index):
        y_index, x_index = divmod(index, self.n_x)
        return x_index, y_index

    def get_next_raw_cell_data(self, delimeter: str = "\t") -> NDArray:
        """
        Load cell data as numpy array using provided cell loader

        :param delimeter: delimeter used in the csv file
        :returns: cell data as numpy array
        """
        return load_cells_from_file_generator(self.files, delimeter)

    def estimate_grid_dimensions(self):
        """
        Estimate how large the row and column spacings need to be for the grid
        """

        for index, cell_data in enumerate(self.get_next_raw_cell_data()):
            mins = np.min(cell_data[:, 2:4], axis=0)
            maxs = np.max(cell_data[:, 2:4], axis=0)
            width, height = maxs - mins

            x_index, y_index = self.get_row_col_from_flat_index(index)
            self.box_widths[x_index] = max(width, self.box_widths[x_index])
            self.box_heights[y_index] = max(height, self.box_heights[y_index])

    def create_tempory_output_files_for_mapped_colonies(self, tmp_dir: str = "./temp"):
        """
        Shift each colony centre by the amount required to fix inside its box.

        Due to the possibility of very large data structures, the temporary states are saved to file
        reduce memory load.
        """
        os.makedirs(tmp_dir, exist_ok=False)
        xs = find_centres(self.box_widths)
        ys = find_centres(self.box_heights)

        for index, cell_data in enumerate(self.get_next_raw_cell_data()):
            x_index, y_index = self.get_row_col_from_flat_index(index)
            cell_data[:, 2:4] += np.array([xs[x_index], ys[y_index]])
            np.savetxt(f"{tmp_dir}/{index:05d}.dat", cell_data)

    def add_all_colonies_to_plot(self, tmp_dir: str = "./temp"):
        """
        Add all the mapped colonies to plot
        """

        files = glob.glob(f"{tmp_dir}/*")
        fig, ax = plt.subplots(1, 1, figsize=ut.getFigsize(cols=2), constrained_layout=True)
        ax.set_facecolor("k")

        for file in files:
            cell_data = np.loadtxt(file)
            nc = len(cell_data)
            cells = [RodShapedBacteria.RodShapedBacterium(ii, *cell_data[ii]) for ii in range(nc)]
            addAllCellsToPlot(cells, ax, ax_rng=1000, alpha=0.8, ec='w')

        fig.savefig(
            f"{tmp_dir}/morphology_phase_diagram.png",
            format="png",
            bbox_inches="tight",
            dpi=400,
        )
        plt.show()


def find_centres(array):
    n = len(array)
    locs = np.zeros(n, dtype=int)
    locs[1:] = array[1:-1] + array[2:]
    return np.cumsum(locs)


if __name__ == "__main__":
    x = np.arange(3)
    y = np.arange(3)
    base_file_name = "/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/SimOutput/ChainingPhaseDiagramFinal/run379/repeat0/"
    files = [f"{base_file_name}/biofilm_{ii:05d}.dat" for ii in range(99, 108)]
    print(len(x), len(y), len(files))
    mpd = MorphologyPhaseDiagram(x, y, files, files)
    mpd.estimate_grid_dimensions()
    try:
        mpd.create_tempory_output_files_for_mapped_colonies()
    except:
        pass
    print(mpd.box_heights)
    print(mpd.box_widths)
    mpd.add_all_colonies_to_plot()
