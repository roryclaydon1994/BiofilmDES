// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <memory>
#include <unordered_map>

// User defined libraries
#include "constants.hpp"     // definition of constants namespace
#include "MathUtility.hpp"
#include "RandUtil.hpp"
#include "Geometry.hpp"
#include "IBacterium.hpp"
#include "SphericalBacteria.hpp"
#include "RodShapedBacteria.hpp"

void SphericalBacterium::divide(std::vector<IBacterium*>& cell_list)
{

  // Proportion of the mother cell volume the first daughter receives
  const double alpha { gen_rand.getNormalRand(0.5,0.1) };
  const double cell1_radius { mRadius*pow(  alpha,1.0/3.0) };
  const double cell2_radius { mRadius*pow(1-alpha,1.0/3.0) };
  const double max_sep { mRadius - ( cell1_radius+cell2_radius ) };

  // Find the RCM of the daughters
  Vec3 rcm_1{
    mPos + Vec3{
      gen_rand.getUniformRand(-0.5*max_sep,0.5*max_sep),
      gen_rand.getUniformRand(-0.5*max_sep,0.5*max_sep),
      gen_rand.getUniformRand(-0.5*max_sep,0.5*max_sep)
    }
  };
  Vec3 rcm_2{
    mPos + Vec3{
      gen_rand.getUniformRand(-0.5*max_sep,0.5*max_sep),
      gen_rand.getUniformRand(-0.5*max_sep,0.5*max_sep),
      gen_rand.getUniformRand(-0.5*max_sep,0.5*max_sep)
    }
  };
#ifndef MOVE_3D // if confined to the plane no randomness in z
  rcm_1.z=0.0;
  rcm_2.z=0.0;
#else
#endif // End move 3d check

  *this = SphericalBacterium
  {
    rcm_1,
    gen_rand.getUniformRand(
      constants::nondim_sphericalGrwthRtePreFac*0.5,
      constants::nondim_sphericalGrwthRtePreFac*1.5
    ),
  cell1_radius
  };
  cell_list.push_back
  (
    new SphericalBacterium
    {
      rcm_2,
      gen_rand.getUniformRand(
        constants::nondim_sphericalGrwthRtePreFac*0.5,
        constants::nondim_sphericalGrwthRtePreFac*1.5
      ),
      cell2_radius
    }
  );
}

#ifdef CHAINING
bool determineLinkedDaughters()
{
 double link_prob { gen_rand.getUniformRand(0,1) };
 // std::cout << "link_prob: " << link_prob << '\n';
 if ( link_prob <= RodShapedBacterium::mLinkingProb ) return true;
 else return false;
}
#endif

void RodShapedBacterium::divide(std::vector<IBacterium*>& cell_list)
{
  // std::cout << "Divide Rod" << '\n';
#if defined(CHAINING)
  // Save the mother links for the daughters to inherit
  IBacterium* mother_upper_link { mUpperEndLinkedTo };
  IBacterium* mother_lower_link { mLowerEndLinkedTo };
#endif

  // Find the rcm of the daughter_cells
  double quarter_full_length{
    0.25*( mLength+2*mRadius )
  };

  Vec3 rcm_1{
    mPos + quarter_full_length * getOrientation()
  };
  Vec3 rcm_2{
    mPos - quarter_full_length * getOrientation()
  };

  *this = RodShapedBacterium
  {
    rcm_1,
    gen_rand.getUniformRand(mAngles.x-1e-3*constants::pi,
                            mAngles.x+1e-3*constants::pi),
#ifndef MOVE_3D // if confined to the plane no randomness in alpha
    0.5*constants::pi,
#else
    gen_rand.getUniformRand(mAngles.y-1e-3*constants::pi,
                            mAngles.y+1e-3*constants::pi),
#endif
    gen_rand.getUniformRand(
      mAvgGrwthRate*0.5,
      mAvgGrwthRate*1.5
    ),
    0.5*(mAvgDivLen-2*mRadius)
  };

  cell_list.push_back(
    new RodShapedBacterium
    {
      rcm_2,
      gen_rand.getUniformRand(mAngles.x-1e-3*constants::pi,
                              mAngles.x+1e-3*constants::pi),
#ifndef MOVE_3D // if confined to the plane no randomness in alpha
      0.5*constants::pi,
#else
      gen_rand.getUniformRand(mAngles.y-1e-3*constants::pi,
                              mAngles.y+1e-3*constants::pi),
#endif
      gen_rand.getUniformRand(
        mAvgGrwthRate*0.5,
        mAvgGrwthRate*1.5
      ),
      0.5*(mAvgDivLen-2*mRadius)
    }
  );

#ifdef CHAINING
  assert( mLowerEndLinkedTo==nullptr );
  assert( mUpperEndLinkedTo==nullptr );
  IBacterium* other_daughter { cell_list.back() };
  assert( other_daughter->getLowerLink()==nullptr );
  assert( other_daughter->getUpperLink()==nullptr );

  // The daughters need to inherit the links from the mother
  mUpperEndLinkedTo=mother_upper_link;
  other_daughter->setLowerLink( mother_lower_link );

  // The cells the daughters are now linked to need to be linked to the daughters
  if ( mother_upper_link ) mother_upper_link->setLowerLink( this );
  if ( mother_lower_link ) mother_lower_link->setUpperLink( other_daughter );

  if ( determineLinkedDaughters() )
  {
    // Create link between the new cells
    mLowerEndLinkedTo = cell_list.back();
    other_daughter->setUpperLink(this);
    // std::cout << "Check cell connections" << '\n';
    // if (mUpperEndLinkedTo)
    // {
    //   std::cout << "A1: " << mUpperEndLinkedTo->getID() << '\n';
    //   std::cout << "B2: " << mUpperEndLinkedTo->getLowerLink()->getID() << '\n';
    // }
    // if (mLowerEndLinkedTo)
    // {
    //   std::cout << "A2: " << mLowerEndLinkedTo->getID() << '\n';
    //   std::cout << "B1: " << mLowerEndLinkedTo->getUpperLink()->getID() << '\n';
    // }
  }
#endif // End chaining

}
