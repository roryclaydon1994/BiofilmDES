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

// std::mt19937 RodShapedBacterium::mGenerator { constants::SEED };
#if defined(ADHESION)
// Set interaction hyperparameters
double RodShapedBacterium::mKappaDep { constants::kappa_depletion };
double RodShapedBacterium::mRd { constants::nondim_Rd };
double RodShapedBacterium::mRi { constants::nondim_Ri };
double RodShapedBacterium::mRepStrength { constants::repulsion_strength };
#elif defined(CHAINING)
double RodShapedBacterium::mKappa { constants::nondim_kappa };   //!< Spring tension
double RodShapedBacterium::mBendRig { constants::nondim_K_bend };//!< Bending rigidity
double RodShapedBacterium::mLinkingProb { 0.5 };                 //!< Probability daughters link
double RodShapedBacterium::mForceThresh { 100 };                 //!< Threshold force before breaking
#elif defined(AG43)
double RodShapedBacterium::mKappa { constants::nondim_kappa };           //!< Spring tension
double RodShapedBacterium::mForceThresh {
  constants::max_extension*constants::nondim_kappa
}; //!< Threshold force before breaking
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

  // const double radius = 0.5*mDiameter;
  // Average area of spherocylinder based on linear growth over lifetime of cell
  // mAvArea = constants::pi*radius*radius + 1.5*radius*mAvgDivLen;

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
{
  // assert(rcm.size()==3);
}

#ifdef ADHESION
void initialiseAdhesionParameters(
  double kappa_depletion,
  double Rd,
  double Ri,
  double rep_strength
)
{
  if ( Rd<constants::nondim_Rd || Rd>=Ri )
  {
    std::cout << "Invalid value of Rd" << '\n';
    exit(10);
  }
  if ( (kappa_depletion<=0) || (Rd<=0) || (Ri<=0) || (rep_strength<=0) )
  {
    std::cout << "Error! All input values must be positive!" << '\n';
    std::cout << "You entered: " << '\n';
    std::cout << "./main.out "
              << kappa_depletion << " "
              << Rd              << " "
              << Ri              << " "
              << rep_strength    << '\n';
    std::cout << "Exiting..." << '\n';
    exit(11);
  }

  RodShapedBacterium::mKappaDep    = kappa_depletion;
  RodShapedBacterium::mRd          = Rd;
  RodShapedBacterium::mRi          = Ri;
  RodShapedBacterium::mRepStrength = rep_strength;
  std::cout << "The hyperparmeters are set as:"   << '\n';
  std::cout << "mKappaDep\t: "
            <<  RodShapedBacterium::mKappaDep     << '\n';
  std::cout << "mRd\t: "
            <<  RodShapedBacterium::mRd           << '\n';
  std::cout << "mRi\t: "
            <<  RodShapedBacterium::mRi           << '\n';
  std::cout << "mRepStr\t: "
            <<  RodShapedBacterium::mRepStrength  << '\n';
}
#endif

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

// RodShapedBacterium::~RodShapedBacterium ()
// {
//   /* Do nothing, vectors are automatically cleared after leaving scope */
// }

// double RodShapedBacterium::consumeNutrients(double c, double c_half)
// {
//   return c / (c + c_half);
// }

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

// void RodShapedBacterium::grow(double c, double dt)
// {
//
//   double area_ratio = getCellArea() / mAvArea;
//   double grwth_rate = mGrwthRtePreFac*area_ratio*consumeNutrients(c);
//   mLength += dt*grwth_rate;
//
// }

// bool RodShapedBacterium::signalDivide()
// {
//   if (mLength>=RodShapedBacterium::mAvgDivLen) {
//     return true;
//   }
//   else return false;
// }

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

// void RodShapedBacterium::setRCMVelFromForce()
// {
//   // exit(EXIT_FAILURE); // How to nondim?
//   // mVel = mForce / (mLength * constants::zeta);
//   mVel = mForce / mLength;
// }

// void RodShapedBacterium::setRCMAngVelFromTorque()
// {
//   // mAngVel = mTorque * 12 / (pow(mLength,3) * constants::zeta);
//   // mAngVel = mTorque * 12 / pow(mLength,3);
//   mAngVel = mTorque * 12 / (mLength*mLength*mLength);
// }

void RodShapedBacterium::reset()
{
  mVel.zero();
  mAngVel.zero();
  mForce.zero();
  mTorque.zero();
}

//
// double RodShapedBacterium::getModE() const
// {
//   // return mRodModE;
//   return constants::nondim_rodModE;
// }

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

// void RodShapedBacterium::setMTSeed(long seed)
// {
//   mSusGenerator.seed(seed);
// }


// double RodShapedBacterium::getDaughterTheta()
// {
//   return getDaughterAngle(mAngles.x);
// }
//
// double RodShapedBacterium::getDaughterAlpha()
// {
//   return getDaughterAngle(mAngles.y);
// }
//
// double getDaughterAngle(double angle)
// {
//   // Assume normally distributed ``kicks'' with small deviation
//   // cell orientations were perturbed with 0.1% uniformly distributed noise
//   // after each division
//   double angle_dist_std_dev{ 1e-3*constants::pi };
//   std::uniform_real_distribution<double> angle_dist(
//     angle-0.5*angle_dist_std_dev,
//     angle+0.5*angle_dist_std_dev
//   );
//   double daughter_angle { angle_dist(RodShapedBacterium::mGenerator) };
//   // std::cout << "angle: " << daughter_angle << '\n';
//   return daughter_angle;
// }
//
// double RodShapedBacterium::getDaughterGrowthRate()
// {
//   std::uniform_real_distribution<double> grwth_dist(
//     constants::nondim_rodGrwthRtePreFac*0.5,
//     constants::nondim_rodGrwthRtePreFac*1.5
//   );
//   double gamma { grwth_dist(RodShapedBacterium::mGenerator) };
//   // std::cout << "gamma: " << gamma << '\n';
//   return gamma;
// }
//
// void createDaughters(
//   uint ii,
//   std::vector<RodShapedBacterium> &cell_list
// )
// {
//   RodShapedBacterium& mother = cell_list[ii];
//
//   // Find the rcm of the daughter_cells
//   double quarter_full_length{
//     0.25*( mother.mLength+2*mother.mRadius )
//   };
//
//   Vec3 rcm_1{
//     mother.mPos + quarter_full_length * mother.getOrientation()
//   };
//   Vec3 rcm_2{
//     mother.mPos - quarter_full_length * mother.getOrientation()
//   };
//
//   mother = RodShapedBacterium {
//     rcm_1,
//     mother.getDaughterTheta(),
//     mother.getDaughterAlpha(),
//     mother.getDaughterGrowthRate()
//   };
//   cell_list.emplace_back(
//     rcm_2,
//     mother.getDaughterTheta(),
//     mother.getDaughterAlpha(),
//     mother.getDaughterGrowthRate()
//   );
//   mother.mAngVel.zero();
//   cell_list.back().mAngVel.zero();
// }

#ifdef AG43
void initialiseAg43Parameters(
  double div_len,
  double grwth_rate,
  double kappa,
  double force_thresh
)
{
  RodShapedBacterium::mAvgDivLen   = div_len;
  RodShapedBacterium::mAvgGrwthRate= grwth_rate;
  RodShapedBacterium::mKappa       = kappa;
  RodShapedBacterium::mForceThresh = force_thresh;

  std::cout << "The hyperparmeters are set as:"   << '\n';
  std::cout << "mAvgDivLen\t: "
            <<  RodShapedBacterium::mAvgDivLen    << '\n';
  std::cout << "mAvgGrwthRate\t:"
            <<  RodShapedBacterium::mAvgGrwthRate << '\n';
  std::cout << "mKappa\t: "
            <<  RodShapedBacterium::mKappa        << '\n';
  std::cout << "mForceThresh\t: "
            <<  RodShapedBacterium::mForceThresh  << '\n';
}

/**
  @brief Attempt to make spring links between bacteria
*/
void createSpringLinks(
  std::vector<IBacterium*> &pars
)
{
  // TODO parallelise
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

      bool stuck_to_this_cell {
        cell->getSprings().contains(neighbour_cell->getID())
      };

      if ( different_cells && stuck_to_this_cell==false )
      {
        // find if in contact
        double s,t;
        Vec3 cv;
        const bool in_contact { checkContact(cell,neighbour_cell,s,t,cv) };

        // if so, create a spring
        if ( in_contact )
        {
          // Make a new spring with the index of the cell to which it points
          // printf("%u -- %u\n",cell->getID(),neighbour_cell->getID());
          cell->getSprings().insert(
            std::make_pair<uint,Springs>(
              neighbour_cell->getID(),
              Springs { cell,neighbour_cell,s,t,cv/cv.norm() }
            )
          );
        }
      }
    }
  }
}

void removeSpringLinks(
  std::vector<IBacterium*> &pars
)
{
  #pragma omp parallel for shared(pars) \
          schedule(static) default(none)
  for ( uint ii=0;ii<pars.size();++ii )
  {
    IBacterium* cell { pars[ii] };
    std::erase_if(
      cell->getSprings(), [](const std::pair<uint,Springs>& c){
        const auto& [key,spring] = c;
        return spring.mRemove;
      }
    );
  }
}

void checkSpringLinks(
  std::vector<IBacterium*> &cells
)
{
  #pragma omp parallel for shared(cells) \
          schedule(static) default(none)
  for ( uint ii=0;ii<cells.size();++ii )
  {
    IBacterium* cell { cells[ii] };
    for ( auto &ss : cell->getSprings() )
    {
      auto& [key,spring] = ss;

      // Get the cell this bacterium is supposed to be connected to
      IBacterium* other_cell { spring.mCellB };
      auto &other_springs { other_cell->getSprings() };

      // check linked cells are the ones expected
      assert( key==other_cell->getID() );
      // for ( auto os : other_springs )
      // {
      //   auto& [key,sp] = os;
      //   std::cout << "check: key: " << key << " "
      //             << sp.mCellA->getID() << " <-> " << sp.mCellB->getID()  << '\n';
      // }
      auto os { other_springs.find(cell->getID()) };
      if ( os==other_springs.end() )
      {
        printf("Error! Partnerless spring\n");
        exit(3);
      }
      assert( os->second.mCellB==spring.mCellA );

      assert( os->second.mT == spring.mS );
      assert( os->second.mS == spring.mT );
      assert( isclose(os->second.mOriLink,-spring.mOriLink) );

    }
  }
}
#endif

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
