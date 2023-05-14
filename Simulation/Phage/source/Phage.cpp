/*******************************************************************************
    \file Phage.cpp
    \brief Details the basic phage class

    @section DESCRIPTION
    BiofilmDES  - a program that simulates a growing colony of microbial cells

    This file contains the implementation of the phage base class.

    Contributing authors:
    Rory Claydon, University of Edinburgh, rory.claydon@ed.ac.uk
    Ricardo Del Rio, University of Edinburgh, s1802059@ed.ac.uk

    \copyright GNU Public License.

    @section LICENSE
    Copyright (2020) The University of Edinburgh.

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
#include <memory>
#include <fstream>
#include <omp.h>

// User defined libraries
#include "constants.hpp"                    // definition of constants namespace
#include "MathUtility.hpp"                  // vector definitions
#include "IBacterium.hpp"                   // Bacteria interface defined here
#include "Phage.hpp"                        // header file
#include "VerletGrid.hpp"

/*----------------------------------------------------------------------------*/
/*!
  \class Phage
  \brief Basic DEM implementation of a point phage particle

  Phage base class defined here for use in agent based rodModElling of biofilms
  and their infections by phage.

  The base class is rodModElled as a point particle which infects bacteria if it
  arrives inside them. Infection events are rodModElled by removing the infected
  susceptible and creating an instance of the infected class.

  Ricardo to do:

    Phage:
      - Write a constructor for this class.
      - Fill in empty functions
      - Other potentially useful functions e.g. overload of ostream
        (c.f. __string__ in python)
      - Will the default destructor be sufficient? Why or why not?

    Testing:
      - Create test script
      - Test diffusion (MSD)
      - Determine safe timestep size/check for interaction along line of motion
        i.e. (phage should not teleport through bacteria)
      - Test infection (create a susceptible particle, then make a function to
        determine if an infection occurs, and if so remove the susceptible and
        and create an infected bacteria)

    Simulations:
      Simulate infection in a static 2D colony
        - Are there param wipe out the bacteria? Phase transition?
        - What is the wave speed of infection?
        - Can we understand the topology of the infection (how is it spreading)?

    Extensions:
      Growing colony
      Chaining colony
      Different phage types
      Mutations in bacteria

  @author Rory Claydon
*/
/*-------------------------- Declare statics ---------------------------------*/

uint Phage::counter { 0 };
double Phage::mDiffusionConst { constants::nondim_phage_diffusion };
double Phage::mLysisPeriod { constants::nondim_lysis_period };
double Phage::mLysisRate { 1.0/constants::nondim_lysis_period };
uint Phage::mBurstSize { constants::burst_size };

/*----------------------------------------------------------------------------*/

Phage::Phage(double rcm_x, double rcm_y, double rcm_z) :
  mPos{rcm_x,rcm_y,rcm_z}
{
  mId = counter++;
}

Phage::Phage(const Vec3 &rcm) : Phage{rcm.x,rcm.y,rcm.z} {}

Phage::~Phage()
{
  /* Do nothing, vectors are automatically cleared after leaving scope */
}

/*--------------------- Update utlities ------------------------------------*/

void infectBacteria(
  Phage* phage,
  IBacterium* host
)
{
#ifdef PHAGE
  phage->setRemove(true);
  if ( host->getMOI()>=1 )
  {
    // Append to the multiplicity of infection (MOI); for now no additional
    // effects from high MOI
    host->updateMOI();
    return;
  }

  // Complete phage start to be produced after the eclipse period
  constexpr double beta {
    std::pow(
      constants::nondim_avg_eclipse_period/constants::nondim_std_eclipse_period,
      2
    )
  };
  const double eclipse_period {
    gen_rand.getGammaRand(beta,constants::nondim_avg_eclipse_period/beta)
  };

  // The lysis period is the total time from infection to burst
  constexpr double avg_phg_generation_period {
    constants::nondim_lysis_period-constants::nondim_avg_eclipse_period
  };

  // phage->getLysisPeriod()
  // Significantly reduced the standard deviation
  constexpr double alpha {
    100*std::pow(
      avg_phg_generation_period/constants::nondim_std_lysis_period,
      2
    )
  };

  const double lysis_period {
    eclipse_period+gen_rand.getGammaRand(alpha,avg_phg_generation_period/alpha)
  };
  host->setLysisTime( lysis_period );

  host->setBurstSize(
    std::ceil(
      constants::phage_production_rate*( lysis_period-eclipse_period )
    )
  );
  // std::cout << "eclipse: " << eclipse_period << '\n';
  // std::cout << "lysis period: " << lysis_period << '\n';
  // std::cout << "eclipse_period<lysis_period: " << (eclipse_period<lysis_period) << '\n';
  // std::cout << "burst: " << host->getBurstSize() << '\n';
  // exit(1);
  assert ( eclipse_period<lysis_period );
  // Use the binomial distribution to set the number of phage released
  // constexpr double prob_phage_produced { 0.5 };
  // host->setBurstSize(
  //   gen_rand.getBinomialRand(
  //     static_cast<int>(std::ceil(phage->getBurstSize()/prob_phage_produced)),
  //     prob_phage_produced
  //   )
  // );

  // Set bacteria to be infected
  host->setInfected();
#else
  std::cout << "phage is not defined! Exiting..." << '\n';
  exit(21);
#endif
}

void infectNeighbours(
  int ii,
  std::vector<IBacterium*> &r_allCells,
  std::vector<Phage*>      &r_allPhage,
  GridCells<IBacterium>    &r_grid
)
{
  // std::cout << "access phage " << ii << '\n';
  Phage* phage { r_allPhage[ii] };
  double max_size
  {
    std::max(
      std::max(phage->getPos().x,phage->getPos().y),
      phage->getPos().z
    )
  };
  // std::cout << "check size " << '\n';
  if ( ((max_size<r_grid.mCellSize.x)
        && (max_size<r_grid.mCellSize.y)
        && (max_size<r_grid.mCellSize.z))==false )
  {
    // if the phage is outside the Verlet grid for cells skip it
    return;
  }
  Uint3 cell_indices { r_grid.getCellIndex(phage->getPos()) } ;
  // loop over all the neighbour particles
  for ( int iz=umax(cell_indices.z-1,0); iz<=umin(cell_indices.z+1,r_grid.mGridSize.z-1); ++iz )
  {
    for ( int iy=umax(cell_indices.y-1,0); iy<=umin(cell_indices.y+1,r_grid.mGridSize.y-1); ++iy )
    {
      for ( int ix=umax(cell_indices.x-1,0); ix<=umin(cell_indices.x+1,r_grid.mGridSize.x-1); ++ix )
      {
        Uint3 neighbour_grid_coords = Uint3(ix,iy,iz);

        // Get the first particle in this node
        // std::cout << "get host" << '\n';
        // if ( r_grid.getParticleCellID(neighbour_grid_coords)>r_grid.mGridCells.size() )
        // {
        //   std::cout << "error, acessing cell outside grid" << '\n';
        //   std::cout << "cell in " << r_grid.getParticleCellID(neighbour_grid_coords) << '\n';
        //   std::cout << "# grid cells" << r_grid.mGridCells.size() << '\n';
        //   std::cout << "neighbour grid coords " << neighbour_grid_coords << '\n';
        //   std::cout << "grid " << r_grid.mCellSize << '\n';
        //   exit(1);
        // }
        IBacterium* potential_host {
          r_grid.mGridCells[r_grid.getParticleCellID(neighbour_grid_coords)].mHead
        };
        // std::cout << "loop over neighbour hosts" << '\n';
        while ( potential_host != nullptr ) {
          // Max separation is head to head
          const double max_sep {
            potential_host->getEffectiveR()
          };
          if (
            dot2( phage->getPos()-potential_host->getPos() ) <= dot2( max_sep )
          )
          {
            // get point to line distance
            // if close enough then infect and remove phage
            // bacteria has a flag that is incremented for the number of
            // infections
            Vec3 lower,upper;
            potential_host->getMyEndVecs(lower,upper);
            double phg_to_host_dist {
              getPointToLineSegDist(phage->getPos(),lower,upper)
            };
            if ( phg_to_host_dist<=potential_host->getRadius() )
            {
              // std::cout << "\ncheck erase: num phage: "
              //           << r_allPhage.size() << '\n';
#pragma omp critical (INFECTION_EVENT)
              {
                infectBacteria(phage,potential_host);
              }
              // std::cout << "check erase: num phage: "
              //           << r_allPhage.size() << '\n';
              // exit(1);
              return;
            }
          }
          // Go down the chain of linked particles in the node
          potential_host = potential_host->getHeadLink();
        }
      }
    }
  }
}

void attemptInfections(
  std::vector<IBacterium*> &r_allCells,
  std::vector<Phage*>      &r_allPhage,
  GridCells<IBacterium>    &r_grid
)
{
#pragma omp parallel for \
        shared(r_allPhage,r_grid,r_allCells) schedule(static) default(none)
  for ( int ii=r_allPhage.size()-1; ii>=0; --ii )
  {
    infectNeighbours(ii,r_allCells,r_allPhage,r_grid);
  }

  // Inherently serial
  for ( int ii=r_allPhage.size()-1; ii>=0; --ii )
  {
    if ( r_allPhage[ii]->signalRemoval() )
    {
      delete r_allPhage[ii];
      // std::cout << "try to delete " << ii << '\n';
      r_allPhage.erase(r_allPhage.begin()+ii);
      // std::cout << "successfully deleted " << ii << '\n';
    }
  }
  // std::cout << "testing parallel phage loop" << '\n';
  // exit(1);
}

struct genPhagePos
{
  Vec3 operator()(double r)
  {
    Vec3 angles {
      gen_rand.getUniformRand(0,2*constants::pi), // phi (azimuthal angle)
#ifdef MOVE_3D
      gen_rand.getUniformRand(0,constants::pi),   // theta (polar angle)
#else
      0.5*constants::pi,                          // theta (polar angle)
#endif
      0                                           // psi (rotation about central axis)
    };
    return gen_rand.getUniformRand()*r*getOrientationFromAngles(angles);
  }

  Vec3 operator()(Vec3 v, double r)
  {
    return gen_rand.getUniformRand()*v+operator()(r);
  }
};

void lyseCell(IBacterium* host, std::vector<Phage*>& phage)
{
#ifdef PHAGE
  Vec3 p1,p2; // The poles of the cell ( not necessarily distinct )
  const Vec3 rcm { host->getPos() };
  host->getMyEndVecs(p1,p2);
  const double r { host->getRadius() }; // Host radius
  const Vec3 v { p2-p1 };    // line segment of the main body
  const Vec3 vv { dot2(v) }; // length squared of the main body
  const uint N { host->getBurstSize() };
  std::vector<Phage*> progeny(N,nullptr);
  if ( vv==0 )
  {
    for ( uint ii=0; ii<N; ++ii ) progeny[ii]=new Phage(rcm+genPhagePos()(r));
  }
  else
  {
    for ( uint ii=0; ii<N; ++ii ) progeny[ii]=new Phage(p1+genPhagePos()(v,r));
  }
  // for ( auto pp : progeny )
  // {
  //   std::cout << "pp: " << pp->getID() << '\n';
  // }
  // exit(1);
  // Add the progeny to the main list
  phage.insert(phage.end(),progeny.begin(),progeny.end());
#else
  std::cout << "phage is not defined! Exiting..." << '\n';
  exit(21);
#endif
}

// void Phage::setVelFromForce(double dt)
// {
//   // TODO: If all phage are the same size, better to make it a static variable or a constant
//   double r_nonD = 1.0;    //TODO set to 0.01d_0 or equivalent , d_0 = bacteria length
//   mVel = mForce / r_nonD;
// }

// void Phage::updatePos2D(double dt)
// {
//   mPosition += dt * mVel;
//   mPosition.mZ = 0;
// }
//
// void Phage::updatePos3D(double dt)
// {
//   this->setVelFromForce(dt);
//   mPosition += dt * mVel;
// }

// bool Phage::infectionEvent(RodShapedBacterium &r_host)
// {
//   // Get two vectors along orientation line passing through CM
//   std::array<Vec3,2> ends {r_host.getMyEndVecs()};
//
//   // Phage infects potential host if it is wihthin the cell
//   double phg_to_host_dist { getPointToLineSegDist(mPos,ends[0],ends[1]) };
//   return ( phg_to_host_dist <= constants::nondim_rodSpheroDiam*0.5 );
// }
//
// void Phage::findThermalForce(double dt)
// {
//   //using https://doi.org/10.1016/j.partic.2013.05.005.
//   //constexpr double A { 1.0 };   //TODO constant to non-dimensionalise, liase with Rory
//   double A {this->mDiffusionConst};
//   std::normal_distribution<double> distribution(0.0, 1.0);
//   const double effective_diff_const { A / sqrt(dt) };
//   Vec3 thermal_force {
//     effective_diff_const*distribution(Particle::mGenerator),
//     effective_diff_const*distribution(Particle::mGenerator),
//     effective_diff_const*distribution(Particle::mGenerator)
//   };
// #ifdef CONFINE_2D
//   thermal_force.mZ = 0;
// #endif // End 2d motion
//
//   //adds thermal force to other forces
//   mForce += thermal_force;
// }
//
// void Phage::getAllForces(double dt) {
//   // This is the same as *this.function(args) without derefrencing
//   this->findThermalForce(dt); //adds brownian force w. normal dist.
// }
//
// // std::ostream &operator<<(std::ostream &r_out, const Phage &r_phage) {
// //   r_out << r_phage.mPosition.mX << '\t'
// //         << r_phage.mPosition.mY << '\t'
// //         << r_phage.mPosition.mZ;
// //   return r_out;
// // }
//
// void printPhageToFile(const std::string file_name,
//                       const Phage::phagevec &r_all_phage) {
//   // std::cout << "Saving to " << file_name << '\n';
//   std::ofstream phage_out{file_name, std::ios::out};
//   if (!phage_out) {
//       std::cerr << "Failed to open file for writing!" << std::endl;
//       exit(EXIT_FAILURE);
//   }
//   // RC TODO: add to static variable
//   std::array<std::string, 3> col_headers{"rcm_x", "rcm_y", "rcm_z"};
//   for (int ii = 0; ii < static_cast<int>(col_headers.size()); ii++) {
//       if (ii == static_cast<int>(col_headers.size()) - 1)
//           phage_out << col_headers[ii];
//       else
//           phage_out << col_headers[ii] << "\t";
//   }
//   phage_out << "\n";
//   for (auto &phage : r_all_phage) {
//       phage_out << *phage << "\n";
//   }
//   phage_out.close();
// }
//
// // void Phage::reset()
// // {
// //   mVel.mX =  mVel.mY = mVel.mZ = 0;
// //   mForce.mX    =  mForce.mY    = mForce.mZ    = 0;
// // }
//
// void calcAllForces(Phage::phagevec &r_allphage, double dt) {
//     //TODO extend implementation to passing phagevec by reference (ask Rory)
//   for ( auto &phage : r_allphage )
//   {
//     phage->getAllForces(dt);
//   }
// }
