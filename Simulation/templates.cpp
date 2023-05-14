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

// Ensure the full Biofilm template definition can be seen
#include "constants.hpp"        // definition of constants namespace
#include "MathUtility.hpp"       // definition of Vec3 class

#include "VerletGrid.hpp"
#include "VerletGrid.cpp"

#include "RodShapedBacteria.hpp"
#include "RodShapedBacteria.cpp"

#include "chaining_susceptible_class.hpp"
#include "chaining_susceptible_class.cpp"

#include "biofilm_class.hpp" // Definitions for biofilm
#include "biofilm_class.cpp" // Don't do this anywhere else

#include "forces.hpp" // Definitions of interactions
#include "forces.cpp"

// Explicit instantiation of the Biofilm<RodShapedBacterium> class
template class Biofilm<RodShapedBacterium>;
template class Biofilm<ChainingRodShapedBacterium>;

template class GridCell<RodShapedBacterium>;
template class GridCell<ChainingRodShapedBacterium>;

// template std::array<Vec3,2> getHertzianTorqueBetweenCellBandCellA(const RodShapedBacterium &cellA,
//                                                                   const RodShapedBacterium &cellB);

// template void printVisulationPairInteraction(const RodShapedBacterium &cellA,
//                                              const RodShapedBacterium &cellB,
//                                              std::string file_name);
