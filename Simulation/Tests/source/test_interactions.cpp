// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <cassert>
#include <random>
#include <chrono>
#include <memory>
#include <numeric>
//#include <execution>

// Custom classes
#include "constants.hpp"         // definition of constants namespace
#include "MathUtility.hpp"
// #include "particle.hpp"
#include "RodShapedBacteria.hpp"
// #include "chaining_susceptible_class.hpp"
#include "infected_class.hpp"
#include "Phage.hpp"
#include "IndexGrid.hpp"
#include "forces.hpp"
#include "HyperParams.hpp"
#include "IO.hpp"

IndexGrid tg;
HyperParams hyperParams;

template < class S >
void bruteForce(std::vector<S> &sus_check)
{
  std::cout << "Brute Force" << '\n';
  const long N { sus_check.size() };
  assert( N>0 );
  for ( long ii=0; ii<N; ++ii )
  {
    Vec3 force  = Vec3(0.0,0.0,0.0);
    Vec3 torque = Vec3(0.0,0.0,0.0);
    for ( long jj=0; jj<N; ++jj )
    {
      if (ii==jj) continue;
      Vec3 tmp_force  = Vec3(0.0,0.0,0.0);
      Vec3 tmp_torque = Vec3(0.0,0.0,0.0);
      pairCollision(
        sus_check[ii],
        sus_check[jj],
        tmp_force,
        tmp_torque
      );
      force+=tmp_force;
      torque+=tmp_torque;
      // if ( ii==131 && dot2(tmp_force)>0 )
      //   std::cout << "131: " << jj << " " << tmp_force << '\n';
    }
    sus_check[ii].mVel = force/sus_check[ii].mLength;
    sus_check[ii].mAngVel = 12*torque/pow(sus_check[ii].mLength,3);
    // if ( ii==131 )
    //   std::cout << "vel : " << sus_check[ii].mId << " " << sus_check[ii].mVel << '\n';
  }
}

template < class S >
void setUp(
  std::vector<S>               &sus,
  std::vector<GenericInfected> &inf,
  std::vector<Phage>    &phg
)
{
  // Read in from file
  populateCellsFromFile(
    "/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/InitialConditions/init_large_1.txt",
    // "/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/SimOutput/main_bacteria_00044.dat",
    sus,
    true
  );
  std::cout << "loaded " << sus.size() << '\n';
  // Set the grid and update the HyperParams
  tg.createCLL(sus,inf,phg);
  return;
}

template < class S >
void testInteractions(std::vector<S> &sus)
{
  std::vector<S> sus_check { sus };
  bruteForce(sus_check);

  std::cout << "CLL Force" << '\n';
  interactParticles(
    sus,
    tg.mCellBegin[IndexGrid::Types::SUS],
    tg.mCellEnd[IndexGrid::Types::SUS]
  );
  for ( int ii=0; ii<sus.size(); ++ii )
  {
    assert( sus_check[ii].mId==sus[ii].mId );
    if ( isclose(sus_check[ii].mVel,sus[ii].mVel)==false )
    {
      std::cout << ii << " checked" << '\n';
      std::cout << "BF: " << sus_check[ii].mId << " " << sus_check[ii].mVel << '\n';
      std::cout << "CL: " << sus[ii].mId << " " << sus[ii].mVel << '\n';
    }
    assert( isclose(sus_check[ii].mVel,sus[ii].mVel) );
  }
}

void testSusInteractions()
{
  std::cout << "Testing Bacteria" << '\n';
  std::vector< RodShapedBacterium >     sus;
  std::vector< GenericInfected > inf;
  std::vector< Phage >    phg;
  setUp(sus,inf,phg);
  testInteractions(sus);
}

// void testChainSusInteractions()
// {
  // std::cout << "Testing Chaining" << '\n';
  // std::vector< GenericInfected > inf;
  // std::vector< Phage >    phg;
  // std::vector< RodShapedBacterium > chaining;
  // setUp(chaining,inf,phg);
  // testInteractions(chaining);
// }

int main(int argc, char const *argv[])
{
  testSusInteractions();
  // testChainSusInteractions();
  return 0;
}
