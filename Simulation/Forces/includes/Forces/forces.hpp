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

#ifndef INTERACTIONS_DEFS
#define INTERACTIONS_DEFS

// Standard libraries
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <array>
#include <vector>
#include <utility>
#include <algorithm>
#include <cassert>
#include <random>
#include <memory>

// User defined
#include "IBacterium.hpp"
#include "RodShapedBacteria.hpp"
#include "constants.hpp"        // definition of constants namespace
#include "MathUtility.hpp"      // definition of Vec3 class
#include "RandUtil.hpp"
#include "IO.hpp"
#include "HyperParams.hpp"

inline void pairCollision( IBacterium *A, const IBacterium *B );
/**<
  \brief Update force and torque on A due to all interactions with B

  @param[in,out] A: Cell to be updated due to interaction
  @param[in] B: Neighbouring cell
  @return Void.

*/

inline double computeAdhesionForceMag(
  const IBacterium *A, const IBacterium *B, double sep
);
/**<
  \brief Find the magnitude of the force due to the presence of EPS in the colony
  @param[in] A: Cell to be updated due to interaction
  @param[in] B: Neighbouring cell
  @param[in] sep: separation between virtual sphere centres
  @param[out] force: magnitude of the force on cell A due to this interaction
  @return Void.
*/

inline double computeHertzianForceMag(
  const IBacterium *A, const IBacterium *B, double sep
);
/**<
  \brief Find the force due to steric repulsion
  @param[in] A: Cell to be updated due to interaction
  @param[in] B: Neighbouring cell
  @param[in] sep: separation between virtual sphere centres
  @param[out] force: force on cell A due to this interaction
  @return Void.
*/

inline void computeHertzianForceTorque(
  const IBacterium *A,
  const IBacterium *B,
  double &sep,
  Vec3 &cA,
  Vec3 &cB,
  Vec3 &force_pos,
  Vec3 &force,
  Vec3 &torque
);
/**<
  \brief Find the force and torque due to steric repulsion
  @param[in] A: Cell to be updated due to interaction
  @param[in] B: Neighbouring cell
  @param[out] sep: separation between virtual sphere centres
  @param[out] force_pos: contact point on the cell surface where the force acts
  @param[out] force/torque: force/torque on cell A due to this interaction
  @return Void.
*/

double getPairHertzianEnergy(
  IBacterium *A,
  const IBacterium *B
);
/**<
  \brief Find the energy due to steric interaction
  @param[in] A: First cell in pair
  @param[in] B: Neighbouring cell
  @return energy.
*/

void getSpringEnergy(
  IBacterium* lower_cell, IBacterium* upper_cell,
  double &bend_energy, double &spring_energy
);
/**<
  \brief Find the energy due to bending spring between a pair
  @param[in] lower_cell: Lower cell in the link
  @param[in] upper_cell: Upper cell in the link
  @param[in] bend_energy: Bending energy due to changes in orientation
  @param[in] spring_energy: Extension energy due to cells pulled apart
  @return Void.
*/

inline void computeAg43Force( IBacterium *cell );
/**<
  @brief For all cells this cell is stuck to, find the sticky force between them.
  @param[in] cell: Cell to be updated due to interaction
  @param[in] neighbours: Cells this cell is stuck to
  @param[out] force/torque: force/torque on cell due to this interaction
  @return Void.
*/

inline void interactSurface( IBacterium *cell );
/**<
  \brief Calculate normal force from hard agar surface and add to cell

  @param[in,out] cell: Cell to be updated due to presence of surface
  @return Void.
*/

inline void interactInterface( IBacterium *cell );
/**<
  \brief Calculate effective "interface" force and add to cell

  Effective gravity to stop cells flying off the agar, although it may lead to
  artificially flat colonies so more elaborate descriptions of the interface are
  possible (but not included here yet).

  @param[in,out] cell: Cell to be updated due to presence of surface
  @return Void.
*/

void polyInteractParticles( std::vector<IBacterium*> &pars );
/**<
  \brief Calculate all forces and torques on cells

  This function uses openMP parallelisation.

  @param[in,out] pars: Cell list updated by all possible interactions
  @return Void.
*/

void computeSingleChainInteraction(
  const IBacterium* lower_cell,
  const IBacterium* upper_cell,
  Vec3 &lower_force,
  Vec3 &upper_force,
  Vec3 &lower_torque,
  Vec3 &upper_torque
);
/**<
  \brief Calculate the force on the lower cell due to the chain link

  @param[in] lower cell: first link in the chain
  @param[in] upper cell: second link in the chain
  @param[out] lower_force,upper_force: force on the respective cell
  @param[out] lower_torque,upper_torque: torque on the respective cell
  @return void
*/


void computeChainingInteractions( IBacterium* cell );
/**<
  \brief Calculate the interactions due to all chains attached to this cell

  @param[in,out] cell: force and torque properties will be updated
  @return void
*/

#endif // end file guard
