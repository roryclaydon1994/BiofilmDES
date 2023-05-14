# Standard libraries
import argparse
# Add parent to path
import os
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, '..')

# Custom modules
from visualise_infected_biofilm import VisInfectedBiofilm
from phage_vis_class import Phage
from infected_vis_class import Infected
from RodShapedBacteria import RodShapedBacterium
from visualise_phage import makeMovie


def calcInfectionProgress(inf_file, inf_start, inf_class=Infected):
    """Calculates both the average distance from origin of all infected and the
     distance of the furthest bacteria from start of infection. (Distance measure = Euclidean)

    Parameters
    ----------
    inf_file : str
        Location of Infected data in time-step.
    inf_start : np.array
        2D positional array where infection started.
    inf_class : custom Infected type
        Class of type Infected for biofilm.
    Returns
    -------
    (float, float)
        A tuple of the max distance of infection and the average distance
        of an infected bacteria from origin of infection at that time-step.

    """
    load_dict = {
        "Infected": [inf_file, inf_class]
    }
    vis_bio = VisInfectedBiofilm(load_dict)
    max_dist = 0
    avg_dist = 0
    if len(vis_bio.species_dict["Infected"]) != 0:
        for cell in vis_bio.species_dict["Infected"]:
            dist = np.linalg.norm(inf_start - cell.rcm[0:2])
            avg_dist += dist
            if dist > max_dist:
                max_dist = dist
        avg_dist *= (1 / len(vis_bio.species_dict["Infected"]))
        return max_dist, avg_dist
    else:
        return 0.0, 0.0


def calcNumPhage(phg_file, phg_class=Phage):
    """Returns the number of phage by the length of the phage species in biofilm"""
    # Set up axes and load only infected as is the only species of interest
    load_dict = {
        "Phage": [phg_file, phg_class]
    }
    vis_bio = VisInfectedBiofilm(load_dict)
    if len(vis_bio.species_dict["Phage"]) != 0:
        return len(vis_bio.species_dict["Phage"])
    else:
        return 0


def analyseWave(base_data_path, upper_index=250):
    """ Creates datafile containing the two distances of interest for the colony"""
    # Start of infection, i.e cm of first infected cell
    infection_start = pd.read_csv(base_data_path + "data/main_infected_00000.dat", delimiter="\t",
                                  index_col="cell_id")
    inf_start_loc = infection_start[["com_vec_x", "com_vec_y"]].iloc[0]

    # Hold infection spread data at each timestep for plotting
    df = pd.DataFrame(index=["infMax", "infAvg"], columns=times)
    for ii in range(upper_index + 1):
        args.inf_file = base_data_path + f"data/main_infected_{ii:05d}.dat"
        # Get two measures of infection at each of 250 time-steps
        current_inf_max, current_inf_avg = calcInfectionProgress(inf_file=args.inf_file, inf_start=inf_start_loc)
        df.loc["infMax"].iloc[ii] = current_inf_max
        df.loc["infAvg"].iloc[ii] = current_inf_avg

    fig, axs = plt.subplots(2, sharex=True, figsize=(10, 8))
    axs[0].plot(times, df.loc["infMax"])
    axs[0].set_title("Distance of furthest away infected bacteria from initial infection site")
    axs[1].plot(times, df.loc["infAvg"])
    axs[1].set_title("Average distance of infected bacteria from initial infection site")
    axs[1].set_xlabel("Time in hours")
    plt.savefig(base_data_path + f"wave_infection/WaveInfection.png", bbox_inches='tight')
    plt.close(fig)
    df.to_csv(base_data_path + f"wave_infection/infection_data.txt", sep="\t")


def calculateSpeeds(base_data_path, upper_index=250):
    """Calculates both measures of speed: -maximum distance from origin of infected bacteria at a time step,
                                          -average distance from origin of infected bacteria at a time step.
       Plots these for each colony at each burst size alongside infection data to spot pattern (ensure
       linear relationship appropiate)"""
    # Assumes linear relationship from start to 10h
    read_data = pd.read_csv(base_data_path + "wave_infection/infection_data.txt", delimiter="\t", index_col=0)
    infMax = read_data.loc["infMax"].astype("float")
    infAvg = read_data.loc["infAvg"].astype("float")

    # Analyse and resample (10 points per lysis period) speeds of infection for that colony
    # V close to just fitting all points, specially for avg infection distance, so optional.
    inf_speeds_max = np.zeros(10)
    inf_speeds_avg = np.zeros(10)
    # Intercepts are needed to plot alongside data to ensure "appropriate" relationship
    intercept_max = 0
    intercept_avg = 0

    if upper_index == 250:
        for i in range(1, 10):
            coeff_max = np.polyfit(times[i:200:10], infMax.iloc[i:200:10], 1)
            coeff_avg = np.polyfit(times[i:200:10], infAvg.iloc[i:200:10], 1)
            inf_speeds_max[i - 1] = coeff_max[0]
            inf_speeds_avg[i - 1] = coeff_avg[0]
            intercept_max += coeff_max[1]
            intercept_avg += coeff_avg[1]
        coeff_max = np.polyfit(times[10:210:10], infMax.iloc[10:210:10], 1)
        coeff_avg = np.polyfit(times[10:210:10], infAvg.iloc[10:210:10], 1)
        inf_speeds_max[9] = coeff_max[0]
        inf_speeds_avg[9] = coeff_avg[0]
        intercept_max += coeff_max[1]
        intercept_avg += coeff_avg[1]
        intercept_max /= 10
        intercept_avg /= 10

        # Final estimates of infection speed
        est_speed_max_inf = np.average(inf_speeds_max)
        est_speed_avg_inf = np.average(inf_speeds_avg)

        fit_max = np.poly1d([est_speed_max_inf, intercept_max])  # The averaged fit to the speed of infection
        fit_avg = np.poly1d([est_speed_avg_inf, intercept_avg])  # The averaged fit to the speed of infection

        fig, axs = plt.subplots(2, sharex=True, figsize=(10, 8))

        axs[0].plot(times[5:200:10], infMax.iloc[5:200:10], "ko", label="original data, 1 point per lysis period")
        axs[0].set_title("Distance of furthest away infected bacteria from initial infection site")
        axs[0].plot(times, fit_max(times),
                    label=f"Linear fit to data (infection speed from max = {est_speed_max_inf:.2f} microns/h)")
        axs[0].legend()

        axs[1].plot(times[5:200:10], infAvg.iloc[5:200:10], "ko", label="original data, 1 point per lysis period")
        axs[1].set_title("Average distance of infected bacteria from initial infection site")
        axs[1].plot(times, fit_avg(times),
                    label=f"Linear fit to data (infection speed from avg = {est_speed_avg_inf:.2f} microns/h)")
        axs[1].set_xlabel("Time in hours")
        axs[1].legend()

        plt.savefig(base_data_path + f"wave_infection/WaveInfection_wfit.png", bbox_inches='tight')
        plt.close(fig)

        return est_speed_max_inf, est_speed_avg_inf

    # Attempt to generalise to any simulation time (larger than 12h)
    else:
        print("Unusual total simulation time")
        # Use only 70% of data, under assumption that edge hasn't been reached
        coeff_max = np.polyfit(times[:round(0.7 * upper_index)], infMax.iloc[:round(0.7 * upper_index)])
        coeff_avg = np.polyfit(times[:round(0.7 * upper_index)], infAvg.iloc[:round(0.7 * upper_index)])
        fit_max = np.poly1d(coeff_max)
        fit_avg = np.poly1d(coeff_avg)

        fig, axs = plt.subplots(2, sharex=True, figsize=(10, 8))

        axs[0].plot(times, infMax)
        axs[0].set_title("Distance of furthest away infected bacteria from initial infection site")
        axs[0].plot(times, fit_max(times),
                    label=f"Linear fit to data (infection speed from max = {coeff_max[0]:.2f} microns/h)")
        axs[0].legend()

        axs[1].plot(times, infAvg)
        axs[1].set_title("Average distance of infected bacteria from initial infection site")
        axs[1].plot(times, fit_avg(times),
                    label=f"Linear fit to data (infection speed from avg = {coeff_avg[0]:.2f} microns/h)")
        axs[1].set_xlabel("Time in hours")
        axs[1].legend()

        plt.savefig(base_data_path + f"wave_infection/WaveInfection_wfit.png", bbox_inches='tight')
        plt.close(fig)

        return coeff_max[0], coeff_avg[0]


def otherAnalysis():
    """Includes any other data analysis in future, e.g number of phage
    NOT CURRENTLY USED AS FORMAT HAS CHANGED"""
    # reflects missing colony 3
    sims = [0, 1, 2, 4, 5, 6, 7, 8, 9, 10]
    for k in sims:
        base_data_path = f"./results/init_large_#.txt-0.0025-250000-10000-5-1-c_lowres/Colony{k}/"

        # Hold infection spread data at each timestep for plotting
        df = pd.DataFrame(index=["phgNum"], columns=times)

        # For now only phage evolution
        for ind, ii in enumerate(range(args.input_index, args.upper_index + 1)):
            args.phg_file = f"{base_data_path}main_phage_{ii:05d}.dat"

            current_phg_num = calcNumPhage(phg_file=args.phg_file)
            df.loc["phgNum"].iloc[ind] = current_phg_num

        plt.plot(times, df.loc["phgNum"], label=f"Colony {k}")

    plt.xlabel("Time in hours")
    plt.ylabel("Number of phage")
    plt.legend()
    plt.show()


def analyse(df, burst_size, base, upper_index=250):
    """ For the first colony in any burst size (or param sweep in future) makes a movie for infection spread.
     Also, oversees the calculation of both speeds (max and avg distance)."""
    colonies = np.arange(1, 11)
    for colony in colonies:
        base_data_path = base + f"/colony|id_{colony}|burst_{burst_size}/"
        if colony == 1:
            if not os.path.exists(base_data_path + "animation/infected_biofilm.mp4"):
                makeMovie(0, upper_index, base_path=base_data_path)
        analyseWave(base_data_path, upper_index)
        speed_max, speed_avg = calculateSpeeds(base_data_path, upper_index)
        dict_of_values = {"burst_size": burst_size, "colony_id": colony, "max_speed": speed_max, "avg_speed": speed_avg}
        df = df.append(dict_of_values, ignore_index=True)
    return df


def findErrors(data, base, burst_size):
    """Calculate errors as SE on mean by finding the 10 colony's calculated infection speeds. Plot alongside average
    of positional data for burst_size across all colonies."""

    # Previously was using average of prior speeds, but Aidan recommended to just use the slope of the
    # averaged data as the speed and still use the same error (as SE of mean of slopes of colonies)
    max_speeds = data["max_speed"]
    estimate_max_speed = np.average(max_speeds)
    estimate_max_speed_err = np.std(max_speeds) / np.sqrt(10)
    print(f"Infection speed according to maximum distance is {estimate_max_speed:.3f} +- {estimate_max_speed_err:.4f}" +
          f" microns/hour for burst size: {burst_size}")
    avg_speeds = data["avg_speed"]
    estimate_avg_speed = np.average(avg_speeds)
    estimate_avg_speed_err = np.std(avg_speeds) / np.sqrt(10)
    print(f"Infection speed according to average distance is {estimate_avg_speed:.3f} +- {estimate_avg_speed_err:.4f}" +
          f" microns/hour for burst size: {burst_size}")

    # hold data for infection waves for 10 colonies at each time point
    df_max = pd.DataFrame(index=[f"infMax_{k}" for k in range(1, 11)], columns=times)
    df_avg = pd.DataFrame(index=[f"infAvg_{k}" for k in range(1, 11)], columns=times)
    for k in range(1, 11):
        data_path = f"{base}/colony|id_{k}|burst_{burst_size}/"
        read_data = pd.read_csv(data_path + "wave_infection/infection_data.txt", delimiter="\t", index_col=0)
        df_max.loc[f"infMax_{k}"][:] = read_data.loc["infMax"]
        df_avg.loc[f"infAvg_{k}"][:] = read_data.loc["infAvg"]

    df_max = df_max.astype("float")
    df_avg = df_avg.astype("float")

    # df.apply works column-wise, i.e averages the locations at each time for the 10 colonies.
    max_inf_dist = df_max.apply(np.average)
    avg_inf_dist = df_avg.apply(np.average)
    err_mean_max = np.std(df_max) / np.sqrt(len(df_max.index))
    err_mean_avg = np.std(df_avg) / np.sqrt(len(df_avg.index))

    fig, axs = plt.subplots(2, sharex=True, figsize=(10, 8))
    axs[0].set_title("Distance of furthest infected from origin")
    axs[0].set_xlabel("Time (hours)")
    axs[0].set_ylabel("Distance (microns)")

    axs[0].errorbar(times[::10], max_inf_dist[::10], fmt="ko", yerr=err_mean_max[::10], capsize=3, markersize=0.3,
                    label="Averaged data (all colonies) with error-bars")

    estimate_max_speed = np.poly1d(
        np.polyfit(times[10:200:10], max_inf_dist.iloc[10:200:10], 1))  # AIDAN  INSPIRED CHANGE

    # Plot line of best fit with slope of speed alongside data
    axs[0].plot(times[10:200:10], estimate_max_speed(times[10:200:10]), "b",
                label=f"Infection speed: {estimate_max_speed.c[0]:.2f} +- {estimate_max_speed_err:.2f} microns/hour")
    axs[0].legend()

    axs[1].set_title("Average distance of infected from origin")
    axs[1].set_xlabel("Time (hours)")
    axs[1].set_ylabel("Distance (microns)")

    axs[1].errorbar(times[::10], avg_inf_dist[::10], fmt="ko", yerr=err_mean_avg[::10], capsize=3, markersize=0.3,
                    label="Averaged data (all colonies) with error-bars")
    # Plot line of best fit with slope of speed alongside data (avg speed slow to get started so intercept corrected
    # to start at the end of first lysis period)

    estimate_avg_speed = np.poly1d(
        np.polyfit(times[10:200:10], avg_inf_dist.iloc[10:200:10], 1))  # AIDAN  INSPIRED CHANGE
    axs[1].plot(times[10:200:10], estimate_avg_speed(times[10:200:10]), "b",
                label=f"Infection speed: {estimate_avg_speed.c[0]:.2f} +- {estimate_avg_speed_err:.2f} microns/hour")
    axs[1].legend()

    plt.savefig(base + "/infection_waves.png", bbox_inches='tight')
    plt.close(fig)

    return estimate_max_speed.c[0], estimate_max_speed_err, estimate_avg_speed.c[0], estimate_avg_speed_err


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog='analyse_infection.py',
        usage='%(prog)s [options] path',
        description='Data analysis script for infection through biofilm',
        epilog='Authors: Ricardo Del Rio, Rory Claydon'
    )
    parser.add_argument("--burst_sizes", required=True,
                        help="burst sizes to look at as a list, e.g 2,3,4,5,10")
    parser.add_argument("-U", "--upper_index", type=int, required=True,
                        help="total number of files to sample")
    parser.add_argument("--dt", type=float, required=False, default=0.0025,
                        help="time-step used in simulation")
    parser.add_argument("--tot_time", type=int, required=False, default=250000,
                        help="total simulation time")
    parser.add_argument("--lysis_period", type=int, required=False, default=10000,
                        help="time taken for cell to lyse")
    parser.add_argument("--diff_const", type=float, required=False, default=1.0,
                        help="diffusion constant in microns²/s")
    args = parser.parse_args()

    # Set parameters used in simulation. NOT USED NOW.
    # Time in zeta/|E| which is in hours, i.e 1 unit -> 5e-5 hours
    # Length in microns
    dt = args.dt
    tot_time = args.tot_time
    lysis_period = args.lysis_period
    diff_const = args.diff_const  # Diffusion const in microns²/s

    num_outputs = tot_time / (0.1 * lysis_period)
    if num_outputs != args.upper_index:
        print(f"WARNING: Might be omitting some data, expected {num_outputs} files and you inputted {args.upper_index}")

    # RC: time interval between plots?
    # RC: I guess you could find the number of lines in the file in the event
    # U is out of bounds?
    # RC: time_interval = tot_time / ( num_outputs-1 ) ?
    time_interval = tot_time / num_outputs
    times = np.arange(0, tot_time + 1, time_interval)
    times *= 5E-5  # Convert time to hours

    burst_sizes = [int(s) for s in args.burst_sizes.split(",")]
    df = pd.DataFrame(columns=["burst_size", "colony_id", "max_speed", "avg_speed"])        # speeds per colony
    df1 = pd.DataFrame(index=burst_sizes, columns=["max_v", "max_v_err", "avg_v", "avg_v_err"]) # speeds per burst_size with errors


    for burst_size in burst_sizes:
        try:
            # Check data exists for that burst size
            base = f"biofilm_infection/burst_{burst_size}"
            if not os.path.isdir(base):
                raise FileNotFoundError

            df = analyse(df, burst_size, base, upper_index=args.upper_index)  # Assume U=250 files outputted
            df1.loc[burst_size] = findErrors(df[df["burst_size"] == burst_size], base, burst_size)

        except FileNotFoundError as e:
            # Make sure wave_infection/ is found alongside data/.
            print("FileNotFoundError: Data doesn't exist, maybe try different burst size")
            raise
            exit(-1)
    df.to_csv("biofilm_infection/infection_speeds.txt", sep="\t")
    df1.to_csv("biofilm_infection/infection_speeds_burst.txt", sep="\t")
