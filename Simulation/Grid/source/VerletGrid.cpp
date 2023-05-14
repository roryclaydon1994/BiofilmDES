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

// Standard libraries
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <memory>
#include <algorithm>

// User defined
// #include "RodShapedBacteria.hpp"
#include "VerletGrid.hpp"

// the side length of the cubes defined by the grid
template < class S >
double GridCell<S>::box_width{0.0};

// Number of elements of each side of the grid defined in VerletGrid.hpp
template < class S >
long GridCell<S>::side_Mx{1};

template < class S >
long GridCell<S>::side_My{1};

template < class S >
long GridCell<S>::side_Mz{1};

// Side lengths of the grid defined in VerletGrid.hpp
template < class S >
double GridCell<S>::length_x{0.0};

template < class S >
double GridCell<S>::length_y{0.0};

template < class S >
double GridCell<S>::length_z{0.0};

/*
  Convert row major format of 3D matrix to 1d index
*/
inline long rowMajorToIndex(long ii, long jj, long kk, long side_Mx,
                                                       long side_My,
                                                       long side_Mz)
{
  return ii + ( jj + kk*side_My ) * side_Mx;
}
