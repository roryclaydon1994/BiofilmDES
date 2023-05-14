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

#ifndef CHAINING_SUSCEPTIBLE_CLASS
#define CHAINING_SUSCEPTIBLE_CLASS

/*
  Chaining RodShapedBacterium particle
*/

// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <random>

// User defined libraries
#include "constants.hpp" // definition of constants namespace
#include "MathUtility.hpp"
#include "RodShapedBacteria.hpp"
#include "HyperParams.hpp"
// #include <Eigen/Dense>

/*----------------------------------------------------------------------------*/
/*
  Bacteria base class defined here for use in agent based rodModElling of biofilms
  The parameter values and simulaton method for the biofilm component of the
  code is strongly based on the following article:
  [1] O. Hallatschek et al., Phys. Rev. Lett., 111, 168101 (2013)
      DOI:https://doi-org.ezproxy.is.ed.ac.uk/10.1103/PhysRevLett.111.168101
  [2] Luca Giomi
  Author: Rory Claydon
*/
/*----------------------------------------------------------------------------*/

// Forward declare existence of Biofilm

class ChainingRodShapedBacterium : public RodShapedBacterium
{
  /*
    Chaining susceptible rodModElled as a spherocylinder with Hertzian forces and
    elastic pole connections after chaining
    All quantities are in non-dimensional units
  */

public:

  uint mUpperEndLinkedTo { BIG };
  uint mLowerEndLinkedTo { BIG };

  // ChainingRodShapedBacterium* mHeadLink { nullptr };
  // std::vector<ChainingRodShapedBacterium*> mMyNeighbourList;

  static std::array<std::string,11> col_headers;

  static double mKappa;
  static double mBendRig;
  static double mLinkingProb;

  ChainingRodShapedBacterium (
    double rcm_x=0, double rcm_y=0, double rcm_z=0,
    double theta = 0,
    double alpha = constants::pi*0.5,
    double grwthPreFac = constants::nondim_rodGrwthRtePreFac,
    double init_length = constants::nondim_init_length,
    uint upper_end_linked_to = BIG,
    uint lower_end_linked_to = BIG
  );

  ChainingRodShapedBacterium (
    const Vec3 &rcm,
    double theta = 0,
    double alpha = constants::pi*0.5,
    double grwthPreFac = constants::nondim_rodGrwthRtePreFac,
    double init_length = constants::nondim_init_length,
    uint upper_end_linked_to = BIG,
    uint lower_end_linked_to = BIG
  );

  bool determineLinkedDaughters() const;

  ~ChainingRodShapedBacterium ();
};

// std::array<Vec3,2> getEndVecs(const ChainingRodShapedBacterium &cell);
/*
  Find the two end point vectors of a ChainingRodShapedBacterium
  Returns:
    vec(P1,P2): the top and bottom end point vectors respectively
*/

std::ostream& operator<< (std::ostream &out, const ChainingRodShapedBacterium &cell);
/*
  Create representation of cell for easy I/O
  This will print in the format
  rcm (rx,ry,rz) orientation (ox,oy,oz) length L \
  top_link_to: cell_ID bottom_link_to: cell_ID
*/

// std::ostream& printCellDataToFile(
//   std::ostream &out,
//   const ChainingRodShapedBacterium *cell
// );
/*
  Fill in later
*/

// void printInteraction(std::ostream &out, const ChainingRodShapedBacterium &cell,
//                       const Vec3 &virt_cen,
//                       const Vec3 &force_on_me,
//                       const Vec3 &torque_on_me);

void getSpringForce(
  const ChainingRodShapedBacterium &lower_cell,
  const ChainingRodShapedBacterium &upper_cell,
  Vec3& force,
  Vec3& torque
);
/*
  Find the chaining force given a cell, which then will try to find if it has an
  upper link. As cells come in pairs, with a top always links to a bottom, only
  need to check if there is an upper cell to which it is linked. Future
  implementations should do this on per chain basis.
  Parameters:

  Effect:

*/

void createDaughters(
  uint ii,
  std::vector<ChainingRodShapedBacterium> &cell_list
);
/*
  Create two daughter cells with slightly perturbed orientation and growth rate
  compared to the mother. The daughters will be chained on either side with the
  mother's original links. Links between daughters are drawn from a random
  distribution.
  The daughters' centres are placed at rcm +/- 0.25*(L+d) with length (L-d)/2
  Effect:
    Emplaces two daughter cells at the end of the passed daughter_cells vector
*/

#endif // End fileguard
