"""
Create a 2D colony morphology
"""

# Standard modules
import os
from typing import Iterable, Optional
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


def load_cells_from_file_generator(
    files: Iterable[str], delimeter: str = "\t"
) -> Iterable[NDArray]:
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

    def __init__(self, x, y, files, labels: Optional[Iterable] = None) -> None:
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
        :param labels: if provided, colour different colonies according to this label and add legend
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

    def get_flat_index_from_row_col(self, x_index: int, y_index: int) -> int:
        """
        Convert the row and column to a flattened index.

        :param x_index: index of the x independent variable
        :param y_index: index of the y independent variable
        :return: flattened index
        """
        return x_index + self.n_x * y_index

    def get_row_col_from_flat_index(self, index: int) -> (int, int):
        """
        Convert the row and column to a flattened index.

        :param index: flattened index
        :return x_index: index of the x independent variable
        :return y_index: index of the y independent variable
        """
        y_index, x_index = divmod(index, self.n_x)
        return x_index, y_index

    def get_next_raw_cell_data(self, delimeter: str = "\t") -> Iterable[NDArray]:
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
        xs = get_midway_steps(self.box_widths)
        ys = get_midway_steps(self.box_heights)

        for index, cell_data in enumerate(self.get_next_raw_cell_data()):
            x_index, y_index = self.get_row_col_from_flat_index(index)
            cell_data[:, 2:4] += np.array([xs[x_index], ys[y_index]])
            np.savetxt(f"{tmp_dir}/{index:05d}.dat", cell_data)

    def add_all_colonies_to_plot(
        self, tmp_dir: str = "./temp", line_width_scale: float = 1000
    ):
        """
        Plot all the mapped colonies on a canvas.

        :param tmp_dir: temporary directory where intermediate files are stored
        :param line_width_scale: modifier for the line width on the bacteria and colony outline
        """

        files = glob.glob(f"{tmp_dir}/*")
        fig, ax = plt.subplots(
            1, 1, figsize=ut.getFigsize(cols=2), constrained_layout=True
        )
        ax.set_facecolor("k")

        for file in files:
            cell_data = np.loadtxt(file)
            nc = len(cell_data)
            cells = (
                RodShapedBacteria.RodShapedBacterium(ii, *cell_data[ii])
                for ii in range(nc)
            )
            addAllCellsToPlot(cells, ax, ax_rng=line_width_scale, alpha=0.8, ec="w")

        ax.axis("scaled")

        fig.savefig(
            f"{tmp_dir}/morphology_phase_diagram.png",
            format="png",
            bbox_inches="tight",
            dpi=400,
        )
        plt.show()


def get_midway_steps(array: NDArray) -> NDArray:
    """
    Given an array of step sizes, determine the position of the centre of each
    step, counting the first to be at 0.

    :param array: array of step sizes
    :returns: positions of each centre step
    """
    locs = np.zeros_like(array)
    locs[1:] = array[1:-1] + array[2:]
    return 0.5 * np.cumsum(locs)


def create_example_morphology_plot(load_directory_path: str):
    """
    Create a morphology phase diagram
    """
    x = np.arange(3)
    y = np.arange(3)
    files = [f"{load_directory_path}/biofilm_{ii:05d}.dat" for ii in range(99, 108)]
    mpd = MorphologyPhaseDiagram(x, y, files)
    mpd.estimate_grid_dimensions()
    try:
        mpd.create_tempory_output_files_for_mapped_colonies()
    except FileExistsError:
        pass
    mpd.add_all_colonies_to_plot()


if __name__ == "__main__":
    generated_output_path = "/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput"
    dir_path = (
        f"{generated_output_path}/SimOutput/ChainingPhaseDiagramFinal/run379/repeat0/"
    )
    create_example_morphology_plot(dir_path)
