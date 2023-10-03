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
#include <iomanip>
#include <cmath>
#include <vector>
#include <cassert>
#include <random>
#include <memory>

// User defined
#include "RodShapedBacteria.hpp"
#include "MathUtility.hpp"
#include "IO.hpp"
#include "constants.hpp"        // definition of constants namespace

uint RodShapedBacterium::counter      { 0 };
double RodShapedBacterium::mRadius    { 0.5*constants::nondim_rodSpheroDiam };
double RodShapedBacterium::mRodModE   { constants::nondim_rodModE };
double RodShapedBacterium::mAvgDivLen { constants::nondim_avg_div_L };
double RodShapedBacterium::mAvgGrwthRate { constants::nondim_rodGrwthRtePreFac };

#if defined(CHAINING)
double RodShapedBacterium::mKappa { constants::nondim_kappa };   //!< Spring tension
double RodShapedBacterium::mBendRig { constants::nondim_K_bend };//!< Bending rigidity
double RodShapedBacterium::mLinkingProb { 0.5 };                 //!< Probability daughters link
double RodShapedBacterium::mForceThresh { 100 };                 //!< Threshold force before breaking
#endif

RodShapedBacterium::RodShapedBacterium (
  double _x,
  double _y,
  double _z,
  double theta,
  double alpha,
  double grwthPreFac,
  double init_length
) :
 mPos{_x,_y,_z}, mAngles {theta,alpha,0}, mLength {init_length},
 mGrwthRtePreFac {grwthPreFac}
{
  mId = counter++;           // increment unique counter
  // Prevent avalanche of divisions
  assert( mLength <= mAvgDivLen );

}

RodShapedBacterium::RodShapedBacterium (const Vec3 &rcm,
                          double theta,
                          double alpha,
                          double grwthPreFac,
                          double init_length) :
             RodShapedBacterium{rcm.x,rcm.y,rcm.z,theta,alpha,
                         grwthPreFac,init_length}
{}

void initialiseRodParameters(
  double aspect_ratio,
  double growth_rate
)
{
  RodShapedBacterium::mAvgDivLen = aspect_ratio;
  RodShapedBacterium::mAvgGrwthRate = growth_rate;
  std::cout << "Rod hyperparmeters are set as:"   << '\n';
  std::cout << "mAvgDivLen\t: "
            <<  RodShapedBacterium::mAvgDivLen    << '\n';
  std::cout << "mAvgGrwthRate\t: "
            << RodShapedBacterium::mAvgGrwthRate  << '\n';
}

#ifdef CHAINING
void initialiseChainingParameters(
  double kappa,
  double bend_rig,
  double linking_prob,
  double force_thresh
)
{
  RodShapedBacterium::mKappa       = kappa;
  RodShapedBacterium::mBendRig     = bend_rig;
  RodShapedBacterium::mLinkingProb = linking_prob;
  RodShapedBacterium::mForceThresh = force_thresh;
  std::cout << "The hyperparmeters are set as:"   << '\n';
  std::cout << "mKappa\t: "
            <<  RodShapedBacterium::mKappa        << '\n';
  std::cout << "mBendRig\t: "
            <<  RodShapedBacterium::mBendRig      << '\n';
  std::cout << "mLinkingProb\t: "
            <<  RodShapedBacterium::mLinkingProb  << '\n';
  std::cout << "mForceThresh\t: "
            <<  RodShapedBacterium::mForceThresh  << '\n';
  std::cout << "mAvgDivLen\t: "
            <<  RodShapedBacterium::mAvgDivLen    << '\n';
  std::cout << "mAvgGrwthRate\t: "
            << RodShapedBacterium::mAvgGrwthRate  << '\n';
}
#endif

double RodShapedBacterium::getCellArea()
{

  double radius { 0.5 };
  // Area of spherocylinder
  return constants::pi*radius*radius + mLength;
}

void RodShapedBacterium::grow(double dt)
{
  mLength += dt*mGrwthRtePreFac;
}

void RodShapedBacterium::move(double dt)
{
  // Get the vector ( dphi / dt, dtheta /dt , dpsi / dt )
  Vec3 bodyFrameAngVel { projectAngVelRod(mAngVel,mAngles) };
#ifndef MOVE_3D
  bodyFrameAngVel.y=0;
  mVel.z=0;
#endif
  // Should phi be mod 2 pi?
  mAngles+=bodyFrameAngVel*dt;
  mPos+=mVel*dt;
}

void RodShapedBacterium::reset()
{
  mVel.zero();
  mAngVel.zero();
  mForce.zero();
  mTorque.zero();
}

Vec3 RodShapedBacterium::getOrientation() const
{
  Vec3 n_hat { getOrientationFromAngles(mAngles) };
#ifndef MOVE_3D
  n_hat.z=0;
#endif
  return n_hat;
}

void RodShapedBacterium::getMyEndVecs(Vec3& p, Vec3& q) const
{
  getEndVecs(*this,p,q);
}

// p is the start, q is the end
void getEndVecs(const RodShapedBacterium &cell, Vec3& p, Vec3& q)
{

  Vec3 cell_dir{cell.getOrientation()*cell.mLength};

  p = cell.mPos - 0.5*cell_dir;
  q = cell.mPos + 0.5*cell_dir;

}

std::ostream& operator<< (std::ostream &out, const RodShapedBacterium &cell)
{
  out << "rcm "       << cell.mPos << " orientation " << cell.getOrientation()
      << " length "   << cell.mLength
      << " diameter " << 2*cell.mRadius
      << " id "       << cell.mId;
  return out;
}

std::ostream& printCellDataToFile(std::ostream &out,
                                  const RodShapedBacterium &cell)
{
  return printBaseCellDataToFile(out,cell);
}
