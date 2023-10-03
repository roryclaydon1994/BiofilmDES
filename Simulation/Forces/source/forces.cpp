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
#include <omp.h>

// User defined
#include "VerletGrid.hpp"
#include "forces.hpp"
#include "RodShapedBacteria.hpp"
#include "constants.hpp"         // definition of constants namespace
#include "MathUtility.hpp"       // definition of Vec3 class

void polyInteractParticles( std::vector<IBacterium*> &pars )
{
  // ------------------ Poly Interactions -----------------
  #pragma omp parallel for shared(pars) \
          schedule(static) default(none)
  for ( uint ii=0;ii<pars.size();++ii )
  {
    IBacterium* cell { pars[ii] };
    for ( auto &neighbour_cell : cell->getNeighbourList() )
    {
      const bool different_cells {
        ( cell->getID()!=neighbour_cell->getID() )
        ||
        ( cell->getMyType()!=neighbour_cell->getMyType() )
      };
      if ( different_cells )
      {
        pairCollision(cell,neighbour_cell);
      }
    }
#ifdef CHAINING
    computeChainingInteractions(cell);
#endif
#ifndef MOVE_3D // If not constrained to 2d add a surface and gravity
    interactSurface(cell);
    interactInterface(cell);
#endif // end move 3d
    cell->setVel();
    cell->setAngVel();
  }
  return;
}

inline double computeHertzianForceMag(
  const IBacterium *A, const IBacterium *B, double sep
)
{
  double sigma { A->getRadius() + B->getRadius() };
  double overlap { std::max(sigma-sep,0.0) }; // Cell overlap
  if ( overlap<=0 ) return 0.0;

  // // Hertzian force magnitude
  double force_mag { overlap*sqrt(overlap) };
  return force_mag;
}

inline void computeHertzianForceTorque(
  const IBacterium *A,
  const IBacterium *B,
  double &sep,
  Vec3 &cA,
  Vec3 &cB,
  Vec3 &force_pos,
  Vec3 &force,
  Vec3 &torque
)
{
  std::cout << "Update force with parallel force" << '\n';
  exit(12);
  getMinDist(A,B,sep,cA,cB);  // Get the minimum distance between the bacteria
  assert(sep>0);

  const Vec3 normal_BA { (cA-cB)/sep }; // Normal direction from B to A

  const double steric_force { computeHertzianForceMag(A,B,sep) };
  force = steric_force*normal_BA;

  // Torque 0 is due to force of B on A
  // lever arm is vector from rod center to virt_sphere center on that rod
  // Take force at the surface of the cell
  force_pos = cA - A->getRadius()*normal_BA;
  torque = cross(force_pos-A->getPos(),force);
}

inline void computeHertzianForceTorque(
  IBacterium *A,
  const IBacterium *B,
  Vec3 &cA,
  Vec3 &cB
)
{
  Vec3 cv = cA-cB;
  double sep=cv.norm();
  assert(sep>0.1);
  Vec3 normal_BA { cv/sep }; // Normal direction from B to A

  double steric_force { computeHertzianForceMag(A,B,sep) };
  Vec3 force = steric_force*normal_BA;

  // Torque 0 is due to force of B on A
  // lever arm is vector from rod center to virt_sphere center on that rod
  // Take force at the surface of the cell
  Vec3 force_pos = cA - A->getRadius()*normal_BA;
  Vec3 torque = cross(force_pos-A->getPos(),force);

  A->addForce(force);
  A->addTorque(torque);
}

inline void interactSurface(IBacterium *cell)
{
  // agar lowest and highest possible values
  constexpr double agar_low {
    constants::nondim_min_z - 0.5*constants::nondim_agar_roughness
  };
  constexpr double agar_high { agar_low+constants::nondim_agar_roughness };

  std::array<Vec3,2> poles;
  cell->getMyEndVecs(poles[0],poles[1]);

  // Reduced quantities for use in the Hertzian force calculation
  const double Rstar { cell->getRadius() };
  const double Estar {
    getEffectiveQ(cell->getModE(),constants::nondim_agar_mod_E)
  };

  const int num_poles { 1 + ( dot2(poles[1]-poles[0])!=0 ) };
  for ( int ii=0; ii<num_poles; ++ii )
  {
    // Height of the agar at the poles indexed by 1,2
    const double local_agar_height {
      gen_rand.getUniformRand(agar_low,agar_high)
    };

    const double sep { poles[ii].z-local_agar_height };

    // Cell pole overlap with agar
    const double overlap { std::max(cell->getRadius()-sep,0.0) };
    if ( overlap<=0 ) return;

    // Note: assumed that cells are essentialy incompressible
    // Rstar is just the cell radius here as the wall is an infinite half plane
    const double force_mag { Estar * overlap * sqrt(Rstar*overlap) };
    const Vec3 force { 0,0,force_mag };

    cell->addForce(force);
    const Vec3 force_pos { poles[ii].x, poles[ii].y, local_agar_height };
    const Vec3 torque = cross(force_pos-cell->getPos(),force);
    cell->addTorque(torque);
  }
}

inline void interactInterface(
  IBacterium *cell
)
{
  const Vec3 gravity { 0,0,-constants::effective_g };
  cell->addForce(gravity);
}

/* === Better parallel handling === */
void pairCollision(IBacterium *A, const IBacterium *B)
{

  double s1,t1,s2,t2; // Line parameter values at which they are closest
  Vec3 c1,c2;         // centres of the closet approach on each line resp.

  // COMs
  Vec3 pos_A { A->getPos() };
  Vec3 pos_B { B->getPos() };

  // find the direction Segment of S1 and S2
  Vec3 v1 = 0.5*A->getLength()*A->getOrientation();
  Vec3 v2 = 0.5*B->getLength()*B->getOrientation();

  // Check line was parallel
  bool par = closestApproachLineSegmentsParallel(pos_A,v1,pos_B,v2,s1,t1,s2,t2);

  c1 = pos_A+v1*s1; // centre of the virtual sphere on cell A
  c2 = pos_B+v2*t1; // centre of the virtual sphere on cell B
  computeHertzianForceTorque(A,B,c1,c2);

  // Parallel segements so two points were returned
  if ( par )
  {
    c1 = pos_A+v1*s2; // centre of the virtual sphere on cell A
    c2 = pos_B+v2*t2; // centre of the virtual sphere on cell B
    computeHertzianForceTorque(A,B,c1,c2);
  }
}

inline double computeHertzianEnergyMag(
  IBacterium *A,
  const IBacterium *B,
  Vec3 cA,
  Vec3 cB
)
{
  Vec3 cv = cA-cB;
  double sep=cv.norm();
  assert(sep>0.1);

  constexpr double sigma { 1.0 };
  double overlap { std::max(sigma-sep,0.0) }; // Cell overlap
  double energy { (2.0/5.0)*overlap*overlap*sqrt(overlap) };
  return energy;
}

/* === Better parallel handling === */
double getPairHertzianEnergy(IBacterium *A, const IBacterium *B)
{

  double energy { 0.0 }; // Energy due to an overlap
  double s1,t1,s2,t2;    // Line parameter values at which they are closest
  Vec3 c1,c2;            // centres of the closet approach on each line resp.

  // COMs
  Vec3 pos_A { A->getPos() };
  Vec3 pos_B { B->getPos() };

  // find the direction Segment of S1 and S2
  Vec3 v1 = 0.5*A->getLength()*A->getOrientation();
  Vec3 v2 = 0.5*B->getLength()*B->getOrientation();

  // Check line was parallel
  bool par = closestApproachLineSegmentsParallel(pos_A,v1,pos_B,v2,s1,t1,s2,t2);

  c1 = pos_A+v1*s1; // centre of the virtual sphere on cell A
  c2 = pos_B+v2*t1; // centre of the virtual sphere on cell B
  energy+=computeHertzianEnergyMag(A,B,c1,c2);

  // Parallel segements so two points were returned
  if ( par )
  {
    c1 = pos_A+v1*s2; // centre of the virtual sphere on cell A
    c2 = pos_B+v2*t2; // centre of the virtual sphere on cell B
    energy+=computeHertzianEnergyMag(A,B,c1,c2);
  }
  return energy;
}

/* === Chaining Bacteria Definitions === */
#ifdef CHAINING
inline Vec3 calcBiNormal(const Vec3 &t0, const Vec3 &t1)
{
  return 2 * cross(t0,t1) / ( t0.norm()*t1.norm() + dot(t0,t1) );
}

inline Vec3 calckb1t0(const Vec3 &t0, const Vec3 &t1)
{
  const Vec3 kb_1 { calcBiNormal(t0,t1) };
  return {
    ( 2 * cross(t1,kb_1) - dot(kb_1,kb_1)*( t1 + t0*(t1.norm()/t0.norm()) ) )
      / ( t0.norm()*t1.norm() + dot(t0,t1) )
  };
}

inline Vec3 calckb1t1(const Vec3 &t0, const Vec3 &t1)
{
  const Vec3 kb_1 { calcBiNormal(t0,t1) };
  return {
    ( -2 * cross(t0,kb_1) - dot(kb_1,kb_1)*( t0 + t1*(t0.norm()/t1.norm()) ) )
      / ( t0.norm()*t1.norm() + dot(t0,t1) )
  };
}

inline void computeSingleChainInteraction(
  const IBacterium* lower_cell,
  const IBacterium* upper_cell,
  Vec3 &lower_force,
  Vec3 &upper_force,
  Vec3 &lower_torque,
  Vec3 &upper_torque
)
{
  assert(upper_cell->getLowerLink()->getID()==lower_cell->getID());
  assert(lower_cell->getUpperLink()->getID()==upper_cell->getID());

  // constexpr double radius { 0.5*constants::nondim_rodSpheroDiam };
  const Vec3 lower_cell_n { lower_cell->getOrientation() };
  const Vec3 upper_cell_n { upper_cell->getOrientation() };

  // The top of the lower cell's head to which the bottom of the spring will attach
  const Vec3 lower_head {
    lower_cell->getPos()
    + lower_cell_n*(
        // radius
        + 0.5*lower_cell->getLength()
      )
  };

  // The tail of the uppers cell's head to which the top of the spring will attach
  const Vec3 upper_tail {
    upper_cell->getPos()
    - upper_cell_n*(
        // radius
        + 0.5*upper_cell->getLength()
      )
  };

  const Vec3 low_to_high { upper_tail-lower_head };
  const double mod_l_to_h { low_to_high.norm() };

  const double inv_mod_l_to_h { 1.0 / mod_l_to_h };
  const Vec3 ux12 { low_to_high * inv_mod_l_to_h };


  /*-------------------*/
  /*   Stiff linking   */
  /*-------------------*/
  const Vec3 lower_rod_base {
    lower_head - lower_cell_n*constants::nondim_rodSpheroDiam
  };
  const Vec3 upper_rod_base {
    upper_tail + upper_cell_n*constants::nondim_rodSpheroDiam
  };

  const Vec3 t0 { lower_head - lower_rod_base };

  // Check there is no performance penalty for this
  const Vec3 t1 { low_to_high };

  const Vec3 t2 { upper_rod_base - upper_tail };

  // define inverse lengths
  const double ds1_inv { 2/( t0.norm() + t1.norm() ) };
  const double ds2_inv { 2/( t1.norm() + t2.norm() ) };

  // Alias for bending rigidity
  const double K { RodShapedBacterium::mBendRig };
  const Vec3 stiff_force0 {
    0.5*K* ( calckb1t0(t0,t1) * ds1_inv )
  };
  const Vec3 stiff_force1 {
    -0.5*K* (
      ( calckb1t0(t0,t1) - calckb1t1(t0,t1) ) * ds1_inv - calckb1t0(t1,t2) * ds2_inv
    )
  };
  const Vec3 stiff_force2 {
    -0.5*K* (
      calckb1t1(t0,t1) * ds1_inv + ( calckb1t0(t1,t2) - calckb1t1(t1,t2)) * ds2_inv
    )
  };
  const Vec3 stiff_force3 {
   -0.5*K* ( calckb1t1(t1,t2) * ds2_inv )
 };

  // Both the forces here are conservative, hence forces on 2 is minus total on 1
  lower_force = stiff_force0 + stiff_force1;
  // upper_force = -lower_force; // Calculate this below after adding spring force

  lower_torque  = cross( lower_rod_base - lower_cell->getPos() , stiff_force0 );
  lower_torque += cross( lower_head     - lower_cell->getPos() , stiff_force1 );
  upper_torque  = cross( upper_tail     - upper_cell->getPos() , stiff_force2 );
  upper_torque += cross( upper_rod_base - upper_cell->getPos() , stiff_force3 );

 //  /*-------------------*/
 //  /*   Spring linking  */
 //  /*-------------------*/
 //
  const Vec3 spring_force {
    RodShapedBacterium::mKappa * ( mod_l_to_h - constants::nondim_rodSpheroDiam ) * ux12
  };

  // Both the forces here are conservative, hence forces on 2 is minus total on 1
  lower_force += spring_force;
  upper_force = -lower_force;

  lower_torque += cross( lower_head - lower_cell->getPos(),  spring_force );
  upper_torque += cross( upper_tail - upper_cell->getPos(), -spring_force );
}

void getSpringEnergy(
  IBacterium* lower_cell, IBacterium* upper_cell,
  double &bend_energy, double &spring_energy
)
{
  // constexpr double radius { 0.5*constants::nondim_rodSpheroDiam };
  const Vec3 lower_cell_n { lower_cell->getOrientation() };
  const Vec3 upper_cell_n { upper_cell->getOrientation() };

  // The top of the lower cell's head to which the bottom of the spring will attach
  const Vec3 lower_head {
    lower_cell->getPos()
    + lower_cell_n*(
        // radius
        + 0.5*lower_cell->getLength()
      )
  };

  // The tail of the uppers cell's head to which the top of the spring will attach
  const Vec3 upper_tail {
    upper_cell->getPos()
    - upper_cell_n*(
        // radius
        + 0.5*upper_cell->getLength()
      )
  };

  const Vec3 low_to_high { upper_tail-lower_head };
  const double mod_l_to_h { low_to_high.norm() };

  const double inv_mod_l_to_h { 1.0 / mod_l_to_h };
  const Vec3 ux12 { low_to_high * inv_mod_l_to_h };


  /*-------------------*/
  /*   Stiff linking   */
  /*-------------------*/
  const Vec3 lower_rod_base {
    lower_head - lower_cell_n*constants::nondim_rodSpheroDiam
  };
  const Vec3 upper_rod_base {
    upper_tail + upper_cell_n*constants::nondim_rodSpheroDiam
  };

  const Vec3 t0 { lower_head - lower_rod_base };

  // Check there is no performance penalty for this
  const Vec3 t1 { low_to_high };

  const Vec3 t2 { upper_rod_base - upper_tail };

  // define inverse lengths
  const double ds1_inv { 2/( t0.norm() + t1.norm() ) };
  const double ds2_inv { 2/( t1.norm() + t2.norm() ) };

  // Alias for bending rigidity
  const double K { RodShapedBacterium::mBendRig };

  bend_energy=0.0;
  bend_energy+=0.5*K*dot2(calcBiNormal(t0,t1))*ds1_inv;
  bend_energy+=0.5*K*dot2(calcBiNormal(t1,t2))*ds2_inv;

  double en_pr_fac { 0.5*RodShapedBacterium::mKappa };
  spring_energy=en_pr_fac*dot2( mod_l_to_h - constants::nondim_rodSpheroDiam );
}

inline
void computeChainingInteractions( IBacterium* cell )
{
  // Temporary variables
  Vec3 lower_force;
  Vec3 lower_torque;

  Vec3 upper_force;
  Vec3 upper_torque;

  if ( cell->getUpperLink() )
  {
    computeSingleChainInteraction(
      cell,
      cell->getUpperLink(),
      lower_force,
      upper_force,
      lower_torque,
      upper_torque
    );
    cell->addForce(lower_force);
    cell->addTorque(lower_torque);
  }

  if ( cell->getLowerLink() )
  {
    computeSingleChainInteraction(
      cell->getLowerLink(),
      cell,
      lower_force,
      upper_force,
      lower_torque,
      upper_torque
    );
    cell->addForce(upper_force);
    cell->addTorque(upper_torque);
  }
}
#endif
