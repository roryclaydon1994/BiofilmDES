- [BiofilmPhageDES](#biofilmphagedes)
  * [Introduction](#introduction)
  * [More code details](#more-code-details)
  * [Documentation](#documentation)
  * [Compilation](#compilation)
  * [Files](#files)
  * [Results](#results)
  * [References](#references)

# BiofilmDES

## Introduction

This repo contains a simulation framework for a discrete element simulation of a bacterial colony and bacteriophage (phage). The susceptible bacteria are represented as rod shaped bacterium such as Escherichia coli and the phage are based on tailed phage such as T4, although in the simulations this are currently point particles for simplicity. Spherical bacteria may also be included.

The biofilm component of the code is strongly based on ref. [1]. The dynamics are assumed to be overdamped due to the low Reynold number environment and interactions between bacterium are based on Hertzian contact forces. The code extends previous work to include novel features, such as elastic links between dividing bacterium. This allows for simulations of chaining bacterial colonies, and also allows us to interpolate between non-chaining and fully chaining colonies. I discuss many aspects of altering the chaining interaction in detail in my thesis, but for general users, the take home is that simple changes to cell-cell interactions can have quite profound influences on the evolution of the colony.

<!-- Work before in [3] used growing, flexible, elastic rods to achieve qualitatively similar patterns seen in chaining Escherichia coli and Bacillus subtilis.  -->


<!-- The second important aspect is to include the presence and effect of phage. Hopefully, spontaneous channel formation will provide an interesting topology for the phage to infect. -->

## Rod shaped interactions

An example of how Hertzian forces between spherocylinders are represented is displayed below. On the right, the method of applying the interaction between rods based on two elastic spheres is demonstrated. On the right, the implementation of parallel interactions is shown. This interaction is likely to be updated in the future!

 <img src="Images/pair_interaction.png" width=500 align=left>
 <img src="Images/line_segment.png" width=1000 align=right>

<div style="page-break-after: always"></div>

## More code details
Euler method with basic updating currently used for simplicity but higher order methods will be implemented later. The collision detection is based on finding the distance betweeen line segments in 3D from [4,5] and computing the Hertzian interaction as in [1] based on this overlap.

At the moment, I use a uniform grid to create linked cell lists (giving O(N) complexity) as due to short range interactions only a few cells need to be checked for interation with any other given cell, although Verlet lists may be more efficient for dense packing (need to check this).

Upon division, there is a chance given by the linking probability that cells will be connected by an elastic rod. This is currently drawn from a uniform distribution. The elastic rod has a bending modulus, $B$, and a compression modulus, K. To account for bending stiffness, there are two elements locked into the poles of the bacteria, moved inside so as not to inadvertently get caught in two bacteria overlapping. The figure below shows the set up.

<p align="center">
 <img src="coupled_ecoli_bending_potential_zoom.png" width=750 align=below>
</p>

## Documentation
See ```main.pdf``` for an explanation of how all the forces are calculated and a bit more background. Apart from the last part on the description of the forces, this is still very much a work in progress.

## Compilation
This package uses CMake (v >= 2.8) for compilation and assumes g++ 9.3.0 or above.
To quickly build and make a binary, call `quickMake.sh` from within the `Simulation/` directory.
This will create a build directory `Simulation/build/`, create a MakeFile with CMake, make the package, and run a simple test to verify if it works.

To do this manually repeat the steps in `quickMake.sh`.
A discrete build directory is highly recommended when using CMake and is encouraged.

The code shouldn't take any command line arguments, but legacy instructions are as follows:
Where the command line inputs are in the order log_file_index, spring_constant and linking_probability.
This is likely to change shortly to include the bending stiffness.  

<!-- ## Files

  ```plotting.nb```:
  Mathematica script for visualising output

  ```output_raw_data.tar.gz```:
  Small example output data for growth of a sitff linked biofilm, essentially 2D. Line 1 column headers are:\
  cell_id	length	diameter	com_vec_x	com_vec_y	com_vec_z	orientation_x	orientation_y	orientation_z *neighbours*	upper_link	lower_link\
  Delimeter is "\t".
  * cell_id: (long) unique identifier for each cell
  * length: (double) current cell length from pole to pole (note this does not include caps! end to end length is current length + diameter)
  * diameter: (double) nondimensionalised cell diameter
  * com_vec_x: (double) center of mass x coordinate (same for y,z suffix)
  * *neighbours*: comma seperated list of cell ids of potential contacts and the potential contacts' potential contacts. This is only outputted if CLOSEST is defined.
  * orientation_x: (double) cell orientation vector x component (same for y,z suffix)
  * upper_link: (long) unique identifier of the cell to which the upper end of the present cell is connected to (None if not connected). Similar for lower_link

```large_sim_high_bending_stiffness.mp4```:
  Example output for large sim (~33000 cells) for high bending stiffness ```K=1```, linking probability ```p=0.995``` and spring constant ```kappa=1.1```.

## Making movies
I've been using ```ffmpeg``` in ```bash``` to make movies by stitching together outputted ```png``` files. These are assumed to be in the ```output``` folder and have the filename pattern ```vis_biofilm_%05d.png```. An example for creating the movie and playing it is included below.

```ffmpeg -r 10 -i output/vis_biofilm_%05d.png -c:v libx264 -vf fps=25 -pix_fmt yuv420p output/vis_biofilm_nondim.mp4 && ffplay output/vis_biofilm_nondim.mp4```

## Results

<p align="center">
 <img src="vis_biofilm_00601_2D.png" width=750 align=below>
</p> -->

## References

1. [Mechanically driven growth of quasi-two dimensional microbial colonies,\
    F.D.C. Farrell, O. Hallatschek, D. Marenduzzo, B. Waclaw,\
    Phys. Rev. Lett. 111, 168101 (2013).](https://doi-org.ezproxy.is.ed.ac.uk/10.1103/PhysRevLett.111.168101)

1.  [Intra-colony channels in E. coli function as a nutrient uptake system. \
     L. M.Rooney, W. B. Amos, P. A. Hoskisson et al.,\
     ISME J 14, 2461–2473 (2020).](https://doi.org/10.1038/s41396-020-0700-9)

1.   [Emergence of active nematics in chaining bacterial biofilms. \
     Y. Yaman, E. Demir, R. Vetter et al.,\
     Nat. Commun. 10, 2285 (2019).](https://doi.org/10.1038/s41467-019-10311-z)

1.  [Three-dimensional distinct element simulation of spherocylinder crystallization.\
    L. Pournin, M. Weber, M. Tsukahara et al.\
    Granul. Matter 7, 119–126 (2005).](https://doi.org/10.1007/s10035-004-0188-4)

1.  [A fast algorithm to evaluate the shortest distance between rods,\
    C. Vega, S. Lago,\
    Comput. Chem., 18(1), 55-59 (1994)](https://doi.org/10.1016/0097-8485(94)80023-5)

1.  I would like to thank Bartlomiej Waclaw from the University of Edinburgh for some very useful discussions on algorithm stability, timestep choice and some    potential optimisations to try out in future
