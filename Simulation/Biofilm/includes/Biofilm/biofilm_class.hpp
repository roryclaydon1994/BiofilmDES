/*******************************************************************************
    BiofilmDES  - a program that simulates a growing colony of microbial cells

    Contributing author:
    Rory Claydon, University of Edinburgh, rory.claydon@ed.ac.uk

    Copyright (2020) The University of Edinburgh.

    The software is based on algorithms described in:

    Mechanically driven growth of quasi-two dimensional microbial colonies,
    F.D.C. Farrell, O. Hallatschek, D. Marenduzzo, B. Waclaw,
    Phys. Rev. Lett. 111, 168101 (2013).

    Three-dimensional distinct element simulation of spherocylinder crystallization.
    Pournin, L., Weber, M., Tsukahara, M. et al.
    Granul. Matter 7, 119–126 (2005).

    A fast algorithm to evaluate the shortest distance between rods,
    C. Vega, S. Lago,
    Comput. Chem., 18(1), 55-59 (1994)

    I would like to thank Bartlomiej Waclaw from Edinburgh University for some
    very useful discussions on algorithm stability, timestep choice and some
    potential optimisations to try out in future.

    This file is part of BiofilmDES.

    BiofilmDES is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BiofilmDES is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    A copy of the GNU General Public License can be found in the file
    License.txt or at <http://www.gnu.org/licenses/>.

    Compilation and run from current directory:
      make && ./biofilm.out 0 1.1 0.95

    Further details in the documentation

*******************************************************************************/

#ifndef BIOFILM_CLASS
#define BIOFILM_CLASS
/*
  Biofilm environment

  Simulates a growing biofilm as a collection of susceptible bacteria, each
  considered to be spherocylinders subject to Hertzian forces. The dynamics
  are overdampled due to the low Reynolds number of the medium.

  Calculating the closest distance between line segments is in the first instance
  based on the algorithm presented in:
  Pournin, L., Weber, M., Tsukahara, M. et al.
  Three-dimensional distinct element simulation of spherocylinder crystallization.
  Granular Matter 7, 119–126 (2005). https://doi.org/10.1007/s10035-004-0188-4

  Carlos Vega, Santiago Lago,
  A fast algorithm to evaluate the shortest distance between rods,
  Computers & Chemistry,
  Volume 18, Issue 1,
  1994,
  Pages 55-59,
  ISSN 0097-8485,
  https://doi.org/10.1016/0097-8485(94)80023-5.
  (http://www.sciencedirect.com/science/article/pii/0097848594800235)

  If necessary I will switch to using:
  David Eberly, Geometric Tools, Redmond WA 98052
  https://www.geometrictools.com/
  and specifically on the description on the Robust Computation Between Line Segments

  Plotting in VMD reference:
  "Humphrey, W., Dalke, A. and Schulten, K., `VMD -Visual Molecular
  Dynamics', J. Molecular Graphics, 1996, vol. 14, pp. 33-38."
*/

// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <memory>

// User defined
#include "VerletGrid.hpp"
#include "RodShapedBacteria.hpp"
#include "chaining_susceptible_class.hpp"
#include "constants.hpp"                    // definition of constants namespace
#include "MathUtility.hpp"

template <class S>
class Biofilm {

public:

  // In future this class needs to be updated to allow for polymorphism.
  // For now, simply use pointers to the only type allowed in the class
  using cellvec  = std::vector< S  >;
  using cellpair = std::array < S,2>;

  cellvec mCellList;                 // vector of susceptibles in biofilm
  IndexGrid mTg;

  Biofilm ();
  /*
    Initiate biofilm with default cell orientation pointing along x-axis.
  */

  Biofilm (cellvec &proto_cells);
  /*
    Initiate a biofilm by supplying a vector of pointers to
    Parameters:
      proto_cells: the cells with which to intiate biofilm
      dt: constant evolution timestep (nondim)
  */

  double size();
  /*
    Calculate the coilony size as in [2]
    Returns:
      sum_i ( cell_length_i + cell_diameter )
  */

  cellvec& getCellsInFilm();
  /*
    Return the vector of susceptibles in this film
  */

  // This should not be a member function
  void getAllForceTorque(
    std::vector< std::unique_ptr<GridCell<S>> > &grid_cells
  );
  /*
    Find the net torque and net force on all members in the colony.
    Assume that the colony is composed of cells which are the same physical
    type i.e. all have same diameters and modulus of elasticity. Further, this
    uses the nondimensional equations for the forces and torques.
    Parameters:
      grid_cells: cells tracking positions of particles for linked cell alg.
    Effect:
      Updates the elements mForce and mTorque for all memebers
      in mCellList
  */

  void getAllHertzianForceTorque(
    std::vector< std::unique_ptr<GridCell<S>> > &grid_cells
  );
  /*
    Find the elastic contact forces between cells. As above assumes that the
    colony is composed of cells which are the same physical type i.e. all have
    same diameters and modulus of elasticity. Further, this uses the
    nondimensional equations for the forces and torques.
    Effect:
      Updates the elements mForce and mTorque for all memebers
      in mCellList only due to contact forces
  */

  void getAllSpringForces();
  /*
    Fill in later.
  */

  void getAllSurfaceForceTorque();
  /*
    Fill in later.
  */

  void getAllGravityForceTorque();
  /*
    Fill in later.
  */

  void handleAllDivisionEvents();
  /*
    Fill in later
  */

  void updateOneTimeStep(double dt);
  /*
    Fill in later
  */

  void evolveBiofilm(double col_size_target = constants::colony_size_target,
                     long seed = constants::SEED);
  /*
    Control the temporal evolution of the biofilm.

    This is the main workhorse of the simulaton, and runs until the biofilm
    reaches a size ( = sum_i (l_i + d) ) > some threshold.
    Parameters:
      col_size_target: if the size of the biofilm exceeds this, stop
      seed: seed for the susceptible rndm number generator
    Effect:
      Run biofilm simulation
  */

  void printBiofilmStateToFile(const std::string filename) const;
  /*
    Fill in later
  */

  void populateCellsFromFile(const std::string file_name);
  /*
      Parameters:
        file_name: file for which to load data from
      Effect:
        populate mCellList to match cells in the file.
  */

  void createLogFile(const std::string file_name,const long seed);
  /*
    Print simulation parameters to file
  */

  ~Biofilm ();
};

template < class S >
void printCellsToFile(const std::string file_name,
                      const typename Biofilm<S>::cellvec &element_list);
/*
  Fill in later
*/

template <class S>
void printVisulationPairInteraction(const S &cellA,
                                    const S &cellB,
                                    std::string="cell_pair_interaction_vis.txt");
/*
  Print to file data all required data for visulatisation of forces and positions
  for two susceptibles. Intended as a debugging function.
  Parameters:
    cellA, cellB: Cells to visualise interation
*/

// template <class S>
// typename Biofilm<S>::cellvec::iterator divideBacteria(
//   Biofilm<S> &bio_film,
//   typename Biofilm<S>::cellvec &daughter_cells,
//   typename Biofilm<S>::cellvec::iterator it
// );
/*
  Fill in later
*/

template <class S>
long preallocCellList(Biofilm<S> &bio_film, double col_size_target);
/*
  Parameters:
    bio_film: the film for which to reserve enough elements
    col_size_target: the target size of the colony
  Effect:
    Increase the capacity of mCellList to hold the number of cells required to
    sum to the target, assuming the cells are of the average size
    i.e. reserve col_size_target / (0.75*div_len + 1)
  Returns:
    expected_element_number: number of elements expected in the final colony
*/

template class Biofilm<RodShapedBacterium>;
template class Biofilm<ChainingRodShapedBacterium>;

#endif // End fileguard
