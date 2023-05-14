#ifndef INFECTED_BIOFILM
#define INFECTED_BIOFILM

// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <memory>

// User defined
#include "VerletGrid.hpp"
#include "RodShapedBacteria.hpp"
#include "infected_class.hpp"
#include "Phage.hpp"
#include "constants.hpp"        // definition of constants namespace
#include "MathUtility.hpp"

class InfectedBiofilm
{

public:

  /*----------------------------- Typedefs -----------------------------------*/
  // TODO: Move to namespace
  using sus_vec  = std::vector< std::unique_ptr<RodShapedBacterium>      >;
  using inf_vec  = std::vector< std::unique_ptr<GenericInfected>  >;
  using phg_vec  = std::vector< std::unique_ptr<Phage>     >;

  /*-------------------------- Species Vectors -------------------------------*/
  sus_vec mSusList;                       //!< vector of susceptibles in biofilm
  inf_vec mInfList;                       //!< vector of infected in biofilm
  phg_vec mPhgList;                       //!< vector of phage in biofilm

  /*-------------------------- Hyperparmeters --------------------------------*/


  /*-------------------------- Constructors ----------------------------------*/

  /**
    \brief Create an empty biofilm
  */
  InfectedBiofilm () {}

  /**
    \brief Fill a biofilm using existing vectors of cells

    @param[in] proto_sus: initial list of bacteria
    @param[in] proto_phg: initial list of phage
    @param[in] proto_inf: intial list of infected
    @param[out] mSusList: assert ownership from proto_sus
    @param[out] mInfList: assert ownership from proto_phg
    @param[out] mPhgList: assert ownership from proto_inf
  */
  InfectedBiofilm (sus_vec &proto_sus, inf_vec &proto_inf, phg_vec &proto_phg) :
    mSusList{ std::move(proto_sus) },
    mInfList{ std::move(proto_inf) },
    mPhgList{ std::move(proto_phg) }
  {}

  ~InfectedBiofilm() {}

  /*------------------------- Evolution --------------------------------------*/
  void updateOneTimeStep(double dt)
  {
    //TODO bring to cpp file.

    /**<
      \relates InfectedBiofilm
      \brief Moves entire static biofilm forward in time by dt. I.e handles
      infection, lysis and movement of phage.

      @param[in]: dt: timestep by which to advance biofilm.
      @param[out]: mSusList: mutated sus_vec with infected cells removed.
      @param[out]: mInfList: mutated inf_vec with infected cells added and
      cells which lyse removed.
      @param[out]: mPhgList: mutated phagevec with new phage produced by lysis
      added and those binding to cells removed.
      @return Void.
    */

    // Kill phage trying to infect already infected (superinfection)
    // static bool debug_flag { false };
    // if ( mPhgList.size()==0 && mInfList.size() == 6 )
    // {
    //   debug_flag = true;
    // }

    // if ( debug_flag ) std::cout << "superinfections" << '\n';
    handleSuperinfectionEvents(mPhgList,mInfList);

    // Infect cells with phage inside (check default empty cellvec allinf works)
    // if ( debug_flag ) std::cout << "infections" << '\n';
    handleAllInfectionEvents(mPhgList,mSusList,mInfList);

    // Check whether any cells are ripe for lysing
    // if ( debug_flag ) std::cout << "lysis" << '\n';
    handleAllLysisEvents(mPhgList, mInfList);

    // Find all the forces on all phage - mutates elements of phage list
    // if ( debug_flag ) std::cout << "forces" << '\n';
    calcAllForces(mPhgList,dt);

    // This could also be the same style as calcAllForces
    for ( auto &phage : mPhgList )
    {
      phage->setVelFromForce(dt);
      phage->move(dt);
      phage->reset();
    }

    for ( auto &infected : mInfList )
    {
      infected->mTimeSinceInfection += dt;
    }
  }

};

#endif
