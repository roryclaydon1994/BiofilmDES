# Standard libraries
import argparse
# Add parent to path
import os
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory as transform

sys.path.insert(0, '..')

# Custom modules
from visualise_infected_biofilm import VisInfectedBiofilm
from RodShapedBacteria import RodShapedBacterium
from find_channels import isPointInCell


def predict_speeds():
    speeds = pd.DataFrame(index=["pred_speed", "max_speed", "min_speed"], columns=["speed", "error"])
    pred_speed = np.zeros(10)
    max_speed = np.zeros(10)
    min_speed = np.zeros(10)

    for i in range(1, 11):

        load_dict = {
            "RodShapedBacterium": [f"initial_conditions/init_large_{i}.txt", RodShapedBacterium]
        }
        vis_bio = VisInfectedBiofilm(load_dict)

        # found colonies have a max radius of around 46 microns^2 so make "continuous" line of points 50 microns long
        # and since diameter is 1 micron then make spacing 0.5 microns to ensure only bacteria which have a negligible
        # overlap with line are "omitted". Also made threshold 0.3 as in worst case there is a 0.6 micron region that
        # the line must fall into (at 0.5microns some points taken as if they belonged to two cells)

        line = np.arange(0, 50.5, 0.5)
        line_x = list(zip(line, np.zeros(101)))
        line_y = list(zip(np.zeros(101), line))
        cells_x = []
        cells_y = []

        # Average expected speed based on x and y directions in colony.
        for point in range(len(line)):
            count = 0
            for index, cell in enumerate(vis_bio.species_dict["RodShapedBacterium"]):
                if isPointInCell(line_x[point], cell, thresh=0.3 * cell.diameter):
                    cells_x.append(cell.rcm)

                if isPointInCell(line_y[point], cell, thresh=0.3 * cell.diameter):
                    cells_y.append(cell.rcm)

        # Remove any appearing multiple lines, i.e overlapping successive points on line
        cells_y = np.unique(np.asarray(cells_y), axis=0)
        cells_x = np.unique(np.asarray(cells_x), axis=0)
        speed_x = max(np.linalg.norm(cells_x, axis=1)) / (len(cells_x) * 0.5)
        speed_y = max(np.linalg.norm(cells_y, axis=1)) / (len(cells_y) * 0.5)

        # Prediction on average speed of infection
        print(
            f"for colony {i} might expect: {(speed_x + speed_y) / 2:.2f} +- {abs(speed_x - speed_y) / 2:.2f} microns/h")
        pred_speed[i - 1] = (speed_x + speed_y) / 2

        avg_length = 0
        for cell in vis_bio.species_dict["RodShapedBacterium"]:
            avg_length += cell.length
        avg_length /= len(vis_bio.species_dict["RodShapedBacterium"])

        lys_period = 0.5  # hours
        # Prediction on max speed (radial orientation)
        print(f"for colony {i} might expect max: {avg_length / lys_period} microns/h")
        max_speed[i - 1] = avg_length / lys_period

        cell_diam = 1
        # Prediction on min speed (concentric rings)
        print(f"for colony {i} might expect max: {cell_diam / lys_period} microns/h")
        min_speed[i - 1] = cell_diam / lys_period

    speeds.loc["pred_speed"]["speed"] = np.mean(pred_speed)
    speeds.loc["pred_speed"]["error"] = np.std(pred_speed) / np.sqrt(len(pred_speed))
    speeds.loc["max_speed"]["speed"] = np.mean(max_speed)
    speeds.loc["max_speed"]["error"] = np.std(max_speed) / np.sqrt(len(max_speed))
    speeds.loc["min_speed"]["speed"] = np.mean(min_speed)
    speeds.loc["min_speed"]["error"] = np.std(min_speed) / np.sqrt(len(min_speed))

    speeds.to_csv("biofilm_infection/speeds_prediction.txt", sep="\t")


def compare_predictions():
    burst_sizes = [2, 3, 4, 5, 10]
    simulation = pd.read_csv("biofilm_infection/infection_speeds_burst.txt", delimiter="\t", index_col=[0])
    prediction = pd.read_csv("biofilm_infection/speeds_prediction.txt", delimiter="\t", index_col=[0])

    fig, ax = plt.subplots(figsize=(15, 10))
    # Predictions
    ax.axhline(y=prediction.loc["pred_speed"][0], color='k', linestyle='--')
    ax.axhline(y=prediction.loc["max_speed"][0], color='k', linestyle='--')
    ax.axhline(y=prediction.loc["min_speed"][0], color='k', linestyle='--')
    ax.errorbar(burst_sizes, simulation["max_v"], yerr=simulation["max_v_err"], capsize=3
                , ecolor="k", label="Simulated speed of furthest infection")
    ax.errorbar(burst_sizes, simulation["avg_v"], yerr=simulation["avg_v_err"], capsize=3
                , ecolor="k", label="Simulated speed of average infection distance")
    ax.text(x=9.55, y=prediction.loc["pred_speed"][0] - 0.11, s='Predicted speed', alpha=0.7, color='k')
    ax.text(x=9.5, y=prediction.loc["max_speed"][0] - 0.11, s='Max possible speed', alpha=0.7, color='k')
    ax.text(x=9.5, y=prediction.loc["min_speed"][0] - 0.11, s='Min. possible speed', alpha=0.7, color='k')
    ax.legend()
    ax.set_title("Speed of bacteriophage infection per burst size")
    ax.set_ylabel("Speed in microns/hour")
    ax.set_xlabel("Burst size (number of phage produced each lysis)")
    custom_ticks = list(ax.get_yticks())
    custom_ticks.append(round(prediction.loc["pred_speed"][0], 2))
    custom_ticks.append(round(prediction.loc["min_speed"][0], 2))
    custom_ticks.append(round(prediction.loc["max_speed"][0], 2))
    ax.set_yticks(custom_ticks)
    plt.show()


if __name__ == "__main__":

    if not os.path.exists("biofilm_infection/speeds_prediction.txt"):
        # Calculate min, max and expected speed of infection from colony make up
        predict_speeds()

    # Compare these predictions with those actually found from simulation of infection
    compare_predictions()
