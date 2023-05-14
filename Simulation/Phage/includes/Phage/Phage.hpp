/*******************************************************************************
    \file Phage.hpp
    \brief Details the basic phage class

    @section DESCRIPTION
    BiofilmDES  - a program that simulates a growing colony of microbial cells

    This file contains the definition of the phage base class.

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

#ifndef PHAGE_CLASS
#define PHAGE_CLASS

// Standard libraries
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <memory>

// User defined libraries
#include "constants.hpp"                    // definition of constants namespace
#include "MathUtility.hpp"                  // vector definitions
#include "IBacterium.hpp"                  // Bacterium interface defined here
// #include "particle.hpp"
//#include "infected_class.hpp"               // Infected defined here, inherits from susceptible
#include "RandUtil.hpp"
#include "VerletGrid.hpp"
#include "Geometry.hpp"

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
/*----------------------------------------------------------------------------*/

class Phage
{

public:

    // enum PhageType { T4 };
    uint mId;                         //!< unique id
    // uint mMyGridCell { constants::INF };

    /*------------------- Set phage hyperparmeters -----------------------------*/
    // uint mType { T4 };
    Vec3 mPos;                        //!< centre of mass vector
    // Vec3 mVel{0,0,0};                 //!< centre of mass velocity
    // Vec3 mForce{0,0,0};               //!< net force experienced by this particle
    // Vec2 mAngles{0.5*constants::pi,0};//!< cell major axis angle to z axis, x axis
    // Vec2 mAngVel{0,0};                //!< angular velocity of the particle in
    // Vec3 mTorque{0,0,0};              //!< net torque experienced by this particle
    // double mLength;                   //!< Length of the phage tail
    // static double mRadius;            //!< particle radius
    static uint   counter;               //!< counts the total number of bacteria
    static double mLysisPeriod;          //!< avg lysis period
    static double mLysisRate;            //!< avg lysis rate
    static uint   mBurstSize;            //!< avg number of progeny per lysis
    static double mDiffusionConst;       //!< diffusion constant
    bool mDead { false };                //!< signals to remove this phage

    bool signalRemoval() const { return mDead; }
    void setRemove(bool remove) { mDead=remove; }
    uint getID() const { return mId; }

    //
    Phage(double rcm_x = 0, double rcm_y = 0, double rcm_z = 0);

    Phage(const Vec3 &rcm);
    /*
    * General purpose constructor
    */

    Vec3 getPos() const { return mPos; }

    virtual ~Phage();

    /*--------------------- Update utlities ------------------------------------*/
    virtual void move(double dt)
    {
      // Diffusion step size
      const double sdev { sqrt(2*dt*mDiffusionConst) };
#ifndef MOVE_3D
      const Vec3 rand_dir{
        gen_rand.getNormalRand(0,sdev),
        gen_rand.getNormalRand(0,sdev),
        0
      };
#else
      const Vec3 rand_dir{
        gen_rand.getNormalRand(0,sdev),
        gen_rand.getNormalRand(0,sdev),
        gen_rand.getNormalRand(0,sdev)
      };
#endif
      mPos+=sdev*rand_dir;
      // Reflecting boundaries for phage at z=-0.5 and 0.5*L
      mPos.z-=std::min(mPos.z+constants::nondim_min_z,0.0);
      mPos.z-=std::max(mPos.z-0.5*constants::max_grid_size,0.0);
      mPos.x=periodicBC(mPos.x,constants::max_grid_size);
      mPos.y=periodicBC(mPos.y,constants::max_grid_size);

      // Ensure the phage are well behaved and stay in their box
      assert(mPos.z>=constants::nondim_min_z);
      assert(mPos.z<=0.5*constants::max_grid_size);
      // if ( std::max(fabs(mPos.x),fabs(mPos.y)) > 0.5*constants::max_grid_size )
      // {
      //   std::cout << std::setprecision(3) << "\nerror, phage outside! "
      //             << mPos << " L: " << (0.5*constants::max_grid_size) << " "
      //             << std::max(fabs(mPos.x),fabs(mPos.y))/(0.5*constants::max_grid_size) << '\n';
      // }
      assert(
        std::max(fabs(mPos.x),fabs(mPos.y)) <= 0.5*constants::max_grid_size
      );
    };

    // ----------------------------------- IO ------------------------------------
    void printToFile(std::ostream &out)
    {
      out << "Phage"           << "\t";
      out << mId               << "\t";
      out << 0                 << "\t";
      out << 0                 << "\t";
      out << mPos              << "\t";
      out << Vec3{0.0,0.0,0.0} << "\n";
    }

    void setVelFromForce(double dt);

    /**<
      \brief Set force experienced by phage due to forces on it

      It is assumed the attribute
      mForce has been updated appropriately e.g. by calling calcAllForces.

      @param[in] dt: simulation timestep
      @param[out] mVelocity: updated with current velocity.
      @return Void.
    */

    // void updatePos2D(double dt);

    /**<
      \brief Simulate phage movement in 2 dimensions.

      It is assumed the attribute
      mForce has been updated appropriately e.g. by calling calcAllForces.

      (Careful with extra factor of 2 and init vel: x = Ft**2/2m + v_0*t)

      @param[in] dt: simulation timestep
      @param[out] mPosition: updated with the new position.
      @return Void.
    */

    // void updatePos3D(double dt);

    /**<
      \brief Simulate phage movement in 3 dimensions.

      It is assumed the attribute
      mForce has been updated appropriately e.g. by calling calcAllForces.

      @param[in] dt: simulation timestep
      @param[out] mPosition: updated with the new position.
      @return Void.
    */

    bool infectionEvent(IBacterium &r_host);

    /**<
      Determine if host will become infected.
      Note for Ricardo: Using pointer to base class so that any derived class from
      base will be usable with this function (albeit interpreted as type base).
      Referred to as polymorphism.

      @param[in] r_host: reference to the potential host
      @return infected_flag: true if the phage is determined to infect this host
    */

    static inline double getLysisPeriod() { return mLysisPeriod; }
    static inline double getLysisRate()   { return mLysisRate;   }
    static inline uint   getBurstSize()   { return mBurstSize;   }

    void findThermalForce(double dt);
    /**<
    \brief Class method version: calc all forces on this particle
    @param[in] dt: simulation timestep
    @param[out] mForce: updated with the brownian force on this particle
    @return Void.
    */

    void getAllForces(double dt);
    /**<
    \brief Class method version: calc all forces on this particle
    @param[in] dt: simulation timestep
    @param[out] mForce force currently experienced by this particle.
    @return Void.
    */

    // void reset();
    // /**<
    // \brief clear all dynamical variables
    // @param[out] mVelocity set to 0
    // @param[out] mForce set to 0
    // @return Void
    // */

    // static std::mt19937& genMT();
    // /**<
    //   \brief return the random number generator for this class
    //   NOTE: needs to return a reference as the mt needs to be advanced in the
    //   caller fn
    //   @returns reference to generator so that it can be updated
    // */
};

void attemptInfections(
  std::vector<IBacterium*> &r_allCells,
  std::vector<Phage*>      &r_allPhage,
  GridCells<IBacterium>    &r_grid
);

void lyseCell(IBacterium* host, std::vector<Phage*>& phage);

void findThermalForce(Phage &r_phage, double dt);
/**<
  \brief Get the total force due to thermal fluctations on this particle
         e.g. from F. Chaumeil and M. Crapper, Particuology, 15, 94-106 (2014)

  Some considerations on DEM for colloids (phage \approx colloids):
  F. Chaumeil and M. Crapper, Particuology, 15, 94-106 (2014)
  DOI: https://doi.org/10.1016/j.partic.2013.05.005

  Note the magnitude of the force depends on the timestep used.

  @param[in]  r_phage: reference to #Phage
  @param[in]  dt: simulation timestep
  @param[out] mForce: updated with the force currently experienced by
              this particle.
  @return Void.

*/


// void calcAllForces(Phage::phagevec &r_allphage, double dt);
/**<
  \relates Phage
  \brief Get the total force on a list of particles e.g. from
         F. Chaumeil and M. Crapper, Particuology, 15, 94-106 (2014)

  Some considerations on DEM for colloids (phage \approx colloids):
  F. Chaumeil and M. Crapper, Particuology, 15, 94-106 (2014)
  DOI: https://doi.org/10.1016/j.partic.2013.05.005

  Main idea for us at the moment will be to implement Brownian motion by
  generating the force first. Note the magnitude of the force depends on the
  timestep used.

  There are other things we can consider including phage-phage interactions
  and Van der Waals for bacteria-phage.

  For now just focus on basic diffusion.

  @param[in] r_allphage: reference to a vector of Phage.
  @param[in] dt: the simulation timestep
  @param[out] mForce: update the mForce attribute of all in the list.
  @returns Void.
*/


// std::ostream &operator<<(std::ostream &r_out, const Phage &r_phage);
/**<
  \relates Phage
  \brief Generate representation of phage for easy I/O

  rx (ry,rz) is the com coordinate in x (y,z)

  @param[in] r_out: outstream to append to
  @param[out] r_out: append out with - rx \t ry \t rz
  @returns r_out
*/

// void printPhageToFile(const std::string file_name,
//                       const Phage::phagevec &r_all_phage);
/**<
  \relates Phage
  \brief Output phage to file

  Note: col_headers currently defined in this function

  @param[in] file_name: name of the output file
  @param[in] r_all_phage: vector of unique pointers to phage in the system
  @returns Void
*/

#endif // end fileguard
