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
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <fstream>

// User defined libraries
#include "constants.hpp" // definition of constants namespace
#include "MathUtility.hpp"
#include "RodShapedBacteria.hpp"
#include "chaining_susceptible_class.hpp"
#include "IO.hpp"
#include "HyperParams.hpp"
// #include <Eigen/Dense>

std::array<std::string,11> ChainingRodShapedBacterium::col_headers
{
  "cell_id\t",
  "length\t",
  "diameter\t",
  "com_vec_x\t",
  "com_vec_y\t",
  "com_vec_z\t",
  "orientation_x\t",
  "orientation_y\t",
  "orientation_z\t",
  "upper_link\t",
  "lower_link"
};

double ChainingRodShapedBacterium::mKappa   { constants::nondim_kappa };
double ChainingRodShapedBacterium::mBendRig { constants::nondim_K_bend };
double ChainingRodShapedBacterium::mLinkingProb { 0.95 };
// double ChainingRodShapedBacterium::mAlpha { constants::nondim_alpha };

ChainingRodShapedBacterium::ChainingRodShapedBacterium (
  double _x, double _y, double _z,
  double theta,
  double alpha,
  double grwthPreFac,
  double init_length,
  uint upper_end_linked_to,
  uint lower_end_linked_to
) : RodShapedBacterium { _x, _y, _z, theta, alpha, grwthPreFac, init_length },
    mUpperEndLinkedTo{ upper_end_linked_to },
    mLowerEndLinkedTo{ lower_end_linked_to }
{
}

ChainingRodShapedBacterium::ChainingRodShapedBacterium (
  const Vec3 &rcm,
  double theta,
  double alpha,
  double grwthPreFac,
  double init_length,
  uint upper_end_linked_to,
  uint lower_end_linked_to
) : RodShapedBacterium { rcm, theta, alpha, grwthPreFac, init_length },
    mUpperEndLinkedTo{ upper_end_linked_to },
    mLowerEndLinkedTo{ lower_end_linked_to }
{
}

ChainingRodShapedBacterium::~ChainingRodShapedBacterium()
{
  // delete mUpperEndLinkedTo;
  // delete mLowerEndLinkedTo;
}

// ChainingRodShapedBacterium* ChainingRodShapedBacterium::getUpperLink() const
// {
//   return mUpperEndLinkedTo;
// }
//
// ChainingRodShapedBacterium* ChainingRodShapedBacterium::getLowerLink() const
// {
//   return mLowerEndLinkedTo;
// }
//
// void ChainingRodShapedBacterium::setUpperLink( ChainingRodShapedBacterium* ptr_linked_cell )
// {
//   mUpperEndLinkedTo = ptr_linked_cell;
// }
//
// void ChainingRodShapedBacterium::setLowerLink( ChainingRodShapedBacterium* ptr_linked_cell )
// {
//   mLowerEndLinkedTo = ptr_linked_cell;
// }

std::ostream& appendID(std::ostream &out, const uint link)
{
  if ( link == BIG )
  {
    out << "None";
  }
  else
  {
    out << link;
  }
  return out;
}

bool ChainingRodShapedBacterium::determineLinkedDaughters() const
{
 std::uniform_real_distribution<double> uniform_disr(0,1);
 double link_prob { uniform_disr(Particle::mGenerator) };
 // std::cout << "link_prob: " << link_prob << '\n';
 if ( link_prob <= ChainingRodShapedBacterium::mLinkingProb ) return true;
 else return false;
}

void createDaughters(
  uint ii,
  std::vector<ChainingRodShapedBacterium> &cell_list
)
{
  ChainingRodShapedBacterium& mother = cell_list[ii];

  // Find the rcm of the daughter_cells
  double quarter_full_length{
    0.25*( mother.mLength+mother.getDiameter() )
  };

  Vec3 rcm_1{
    mother.mPos + quarter_full_length * mother.getOrientation()
  };
  Vec3 rcm_2{
    mother.mPos - quarter_full_length * mother.getOrientation()
  };

  uint u1 { mother.mUpperEndLinkedTo };
  uint u2 { BIG };
  uint l1 { BIG };
  uint l2 { mother.mLowerEndLinkedTo };
  if ( mother.determineLinkedDaughters() )
  {
    l1=cell_list.size();  // Original cell links to the new one at the end of the list
    u2=ii;                // New cell links to the original
  }

  mother = ChainingRodShapedBacterium {
    rcm_1,
    mother.getDaughterTheta(),
    mother.getDaughterAlpha(),
    mother.getDaughterGrowthRate(),
    constants::nondim_init_length,
    u1,
    l1
  };
  cell_list.emplace_back(
    rcm_2,
    mother.getDaughterTheta(),
    mother.getDaughterAlpha(),
    mother.getDaughterGrowthRate(),
    constants::nondim_init_length,
    u2,
    l2
  );
  mother.mId=ii;
  cell_list.back().mId=cell_list.size()-1;

  // // Go to the cells in the linked list and update the pointers so that they
  // // now appropriately point to the new daughter cells and not the deceased mother
  //
  // // The mother's upper link needs to set it's lower link to the first daughter
  if ( u1 != BIG )
    cell_list[u1].mLowerEndLinkedTo = u2;
  //
  // // The mother's lower link needs to set it's upper link to the second daughter
  if ( l2 != BIG )
    cell_list[l2].mUpperEndLinkedTo = l1;

}

std::ostream& operator<< (std::ostream &out, const ChainingRodShapedBacterium &cell)
{
  out << static_cast<RodShapedBacterium>(cell) << " ";
  out << "upper_link_ID: "; appendID( out, cell.mUpperEndLinkedTo ) << " ";
  out << "lower_link_ID: "; appendID( out, cell.mLowerEndLinkedTo );
  return out;
}

std::ostream& printCellDataToFile(std::ostream &out,
                                  const ChainingRodShapedBacterium &cell)
{
  printBaseCellDataToFile(out,cell)        << "\t";
  appendID( out, cell.mUpperEndLinkedTo )  << "\t";
  appendID( out, cell.mLowerEndLinkedTo );
  return out;
}

inline Vec3 calcBiNormal(const Vec3 &t0, const Vec3 &t1)
{
  return 2 * t0.cross(t1) / ( t0.norm()*t1.norm() + t0.dot(t1) );
}

inline Vec3 calckb1t0(const Vec3 &t0, const Vec3 &t1)
{
  const Vec3 kb_1 { calcBiNormal(t0,t1) };
  return {
    ( 2 * t1.cross(kb_1) - kb_1.dot(kb_1)*( t1 + t0*(t1.norm()/t0.norm()) ) )
      / ( t0.norm()*t1.norm() + t0.dot(t1) )
  };
}

inline Vec3 calckb1t1(const Vec3 &t0, const Vec3 &t1)
{
  const Vec3 kb_1 { calcBiNormal(t0,t1) };
  return {
    ( -2 * t0.cross(kb_1) - kb_1.dot(kb_1)*( t0 + t1*(t0.norm()/t1.norm()) ) )
      / ( t0.norm()*t1.norm() + t0.dot(t1) )
  };
}

// Finds the force on the first particle
// inline Vec3 calcTermCurvForce(const Vec3 &t0, const Vec3 &t1, const Vec3 &t2)
// {
  // // Avoid taking sqrts more than necessary
  // const double mod_t0 { t0.norm() };
  // const double mod_t1 { t1.norm() };
  // const double mod_t2 { t2.norm() };
  //
  // // Define these denominators which come up a few times
  // const double denom_01 { 1/( mod_t0*mod_t1 + t0.dot(t1) ) };
  // const double denom_12 { 1/( mod_t1*mod_t2 + t1.dot(t2) ) };
  //
  // // Find the curvature binormal between these vectors
  // const Vec3 kb_1 {
  //   2 * t0.cross(t1) * denom_01
  // };
  // const Vec3 kb_2 {
  //   2 * t1.cross(t2) * denom_12
  // };
  //
  // // std::cout << "K: "<< constants::nondim_K_bend << '\n';
  //
  // // Find the forces
  // const Vec3 force_1 {
  //   ( 2*constants::nondim_K_bend / (  mod_t0 + mod_t1 ) )*(
  //     2*t1.cross(kb_1) - kb_1.dot(kb_1) * ( ( t1 + (mod_t1/mod_t0)*t0 ) )
  //   ) * denom_01
  // };
  // std::cout << "force 1" << force_1 << '\n';
  //
  // const Vec3 force_2 {
  //   ( 2*constants::nondim_K_bend / ( mod_t1 + mod_t2 ) )*(
  //     -2*kb_2.cross(t2) - kb_2.dot(kb_2) * ( ( t2 + (mod_t2/mod_t1)*t1 ) )
  //   ) * denom_12
  // };
  // std::cout << "force 2" << force_2 << '\n';
  // const Vec3 force_1 {
  //   ( 2*constants::nondim_K_bend / (  mod_t0 + mod_t1 ) )*(
  //     2*kb_1.cross(t0) - kb_1.dot(kb_1) * t0
  //   ) * denom_01
  // };
  //
  // const Vec3 force_2 {
  //   ( 2*constants::nondim_K_bend / ( mod_t1 + mod_t2 ) )*(
  //     -2*kb_2.cross(t2) - kb_2.dot(kb_2) * t2
  //   ) * denom_12
  // };

//   return force_1 + force_2;
// }

/*
  Return the force and torque on the lower cell
*/
void getSpringForce(
  const ChainingRodShapedBacterium &lower_cell,
  const ChainingRodShapedBacterium &upper_cell,
  Vec3& force,
  Vec3& torque
)
{

  // constexpr double radius { 0.5*constants::nondim_rodSpheroDiam };
  const Vec3 lower_cell_n { lower_cell.getOrientation() };
  const Vec3 upper_cell_n { upper_cell.getOrientation() };

  // The top of the lower cell's head to which the bottom of the spring will attach
  const Vec3 lower_head {
    lower_cell.mPos
    + lower_cell_n*(
        // radius
        + 0.5*lower_cell.mLength
      )
  };

  // The tail of the uppers cell's head to which the top of the spring will attach
  const Vec3 upper_tail {
    upper_cell.mPos
    - upper_cell_n*(
        // radius
        + 0.5*upper_cell.mLength
      )
  };

  const Vec3 low_to_high { upper_tail-lower_head };
  const double mod_l_to_h { low_to_high.norm() };

  const double inv_mod_l_to_h { 1.0 / mod_l_to_h };
  const Vec3 ux12 { low_to_high * inv_mod_l_to_h };

  /*-------------------*/
  /*   Stiff linking   */
  /*-------------------*/
 //  const Vec3 lower_rod_base {
 //    lower_head - lower_cell_n*constants::nondim_rodSpheroDiam
 //  };
 //  const Vec3 upper_rod_base {
 //    upper_tail + upper_cell_n*constants::nondim_rodSpheroDiam
 //  };
 //
 //  const Vec3 t0 { lower_head - lower_rod_base };
 //
 //  // Check there is no performance penalty for this
 //  const Vec3 t1 { low_to_high };
 //
 //  const Vec3 t2 { upper_rod_base - upper_tail };
 //
 //  // define inverse lengths
 //  const double ds1_inv { 2/( t0.norm() + t1.norm() ) };
 //  const double ds2_inv { 2/( t1.norm() + t2.norm() ) };
 //
 //  // Alias for bending rigidity
 //  const double K { ChainingRodShapedBacterium::mBendRig };
 //  const Vec3 stiff_force0 {
 //    0.5*K* ( calckb1t0(t0,t1) * ds1_inv )
 //  };
 //  const Vec3 stiff_force1 {
 //    -0.5*K* (
 //      ( calckb1t0(t0,t1) - calckb1t1(t0,t1) ) * ds1_inv - calckb1t0(t1,t2) * ds2_inv
 //    )
 //  };
 //  const Vec3 stiff_force2 {
 //    -0.5*K* (
 //      calckb1t1(t0,t1) * ds1_inv + ( calckb1t0(t1,t2) - calckb1t1(t1,t2)) * ds2_inv
 //    )
 //  };
 //  const Vec3 stiff_force3 {
 //   -0.5*K* ( calckb1t1(t1,t2) * ds2_inv )
 // };
 //
 //  // Both the forces here are conservative, hence forces on 2 is minus total on 1
 //  force = stiff_force0 + stiff_force1;
 //  // force_pair[1] = -force_pair[0];
 //
 //  torque  = ( lower_rod_base - lower_cell.mPos ).cross( stiff_force0 );
 //  torque += ( lower_head     - lower_cell.mPos ).cross( stiff_force1 );
 //  // torque_pair[1]  = ( upper_tail     - upper_cell.mPos ).cross( stiff_force2 );
 //  // torque_pair[1] += ( upper_rod_base - upper_cell.mPos ).cross( stiff_force3 );

  /*-------------------*/
  /*   Spring linking  */
  /*-------------------*/

  const Vec3 spring_force {
    ChainingRodShapedBacterium::mKappa * ( mod_l_to_h - constants::nondim_rodSpheroDiam ) * ux12
  };

  // Both the forces here are conservative, hence forces on 2 is minus total on 1
  force = spring_force;
  std::cout << "spring force: " << spring_force << '\n';
  // force_pair[1] = -force_pair[0];

  torque = ( lower_head - lower_cell.mPos ).cross(  spring_force );
  std::cout << "spring torque: " << spring_force << '\n';
  // torque_pair[1] += ( upper_tail - upper_cell.mPos ).cross( -spring_force );

}


// template <>
// void Biofilm<ChainingRodShapedBacterium>::getAllSpringForces()
// {
//   /*
//     The way this is currently set up, a top can only link to a bottom so we
//     only need to loop through and calculate all the forces of a top to a bottom.
//   */
//   for ( auto &cell : mCellList )
//   {
//     // If there is no upper link this is a chain head.
//     if ( cell->mUpperEndLinkedTo == nullptr)
//     {
//       continue;
//     }
//
//     std::array<Vec3,2> torque_pair;
//     std::array<Vec3,2> force_pair;
//
//     getSpringForce(cell, cell->mUpperEndLinkedTo, force_pair, torque_pair);
//
//     cell->mForce += force_pair[0];
//     cell->mTorque += torque_pair[0];
//
//     (cell->mUpperEndLinkedTo)->mForce += force_pair[1];
//     (cell->mUpperEndLinkedTo)->mTorque += torque_pair[1];
//
//     // std::cout << "spring force"  << force_pair[0]  << " at " << mCellList.size()<< '\n';
//     // std::cout << "spring torque" << torque_pair[0] << " at " << mCellList.size()<< '\n';
//     // if ( abs(force_pair[0].mZ) > 1e-2 )
//     // {
//     //   Vec3 orientation { cell->getOrientation() };
//     //   std::cout << orientation << " norm: " << orientation.norm() << '\n';
//     //   std::cout << *cell << '\n';
//     //   std::cout << *cell->mUpperEndLinkedTo << '\n';
//     //   exit(EXIT_FAILURE);
//     // }
//   }
//
// }
