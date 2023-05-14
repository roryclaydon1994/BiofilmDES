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
// //#include <execution>

// Custom classes
#include "constants.hpp"         // definition of constants namespace
#include "MathUtility.hpp"
// #include "particle.hpp"
#include "IBacterium.hpp"
#include "SphericalBacteria.hpp"
#include "RodShapedBacteria.hpp"
// #include "chaining_susceptible_class.hpp"
// #include "infected_class.hpp"
// #include "Phage.hpp"
#include "VerletGrid.hpp"
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
  const ulong N { sus_check.size() };
  assert( N>0 );
  for ( ulong ii=0; ii<N; ++ii )
  {
    Vec3 force  = Vec3(0.0,0.0,0.0);
    Vec3 torque = Vec3(0.0,0.0,0.0);
    for ( ulong jj=0; jj<N; ++jj )
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
    }
    sus_check[ii].mForce = force;
    sus_check[ii].mVel = force/sus_check[ii].mLength;
    sus_check[ii].mTorque = torque;
    sus_check[ii].mAngVel = 12*torque/pow(sus_check[ii].mLength,3);
  }
}

template < class S >
void setUp(
  std::vector<S> &sus
  // std::vector<GenericInfected> &inf,
  // std::vector<Phage>    &phg
)
{
  // Read in from file
  populateCellsFromFile(
    "/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/InitialConditions/init_large_1.txt",
    // "/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/SimOutput/main_bacteria_00044.dat",
    sus,
    true
  );
  // Set the grid and update the HyperParams

  return;
}

void cLL(std::vector<IBacterium*> &poly)
{
  // TODO: this needs to be converted to an array of pointers
  // TODO: interactions should be parallelised
  GridCells grid{ poly, constants::box_width };
  std::cout << "mGridSize "     << grid.mGridSize     << "\n"
            << "mNumGridCells " << grid.mNumGridCells << "\n"
            << "mCellSize "     << grid.mCellSize     << "\n"
            << "mOrigin "       << grid.mOrigin       << "\n";
  grid.updateGridCellLists(poly);
  grid.createNeighbourLists(poly);
  std::cout << "CLL Force" << '\n';
  polyInteractParticles(
    poly
  );
}

template < class S >
void testInteractions(std::vector<S> &sus)
{
  std::vector<S> sus_check { sus };
  bruteForce(sus_check);

  std::vector<IBacterium*> poly(sus.size());
  for ( uint ii=0; ii<sus.size(); ++ii )
  {
    poly[ii] = &sus[ii];
  }
  std::cout << "loaded: " << poly.size() << '\n';
  cLL(poly);

  for ( uint ii=0; ii<sus.size(); ++ii )
  {
    assert( sus_check[ii].mId==sus[ii].mId );
    if ( isclose(sus_check[ii].mVel,sus[ii].mVel)==false )
    {
      std::cout << ii << " checked" << '\n';
      std::cout << "BF: " << sus_check[ii] << " " << sus_check[ii].mVel << '\n';
      std::cout << "CL: " << sus[ii] << " " << sus[ii].mVel << '\n';
    }
    assert( isclose(sus_check[ii].mVel,sus[ii].mVel) );
    assert( isclose(sus_check[ii].mForce,sus[ii].mForce) );
    if ( isclose(sus_check[ii].mAngVel,sus[ii].mAngVel)==false )
    {
      std::cout << ii << "/" << sus.size() << " checked" << '\n';
      std::cout << std::setprecision(8)
                << "BF: " << sus_check[ii].mId << " angvel: "
                << sus_check[ii].mAngVel
                << " torque: "
                <<  sus_check[ii].mTorque
                << '\n';
      std::cout << std::setprecision(8)
                << "CL: " << sus[ii].mId  << " angvel: "
                << sus[ii].mAngVel
                << " torque: "
                <<  sus[ii].mTorque
                << '\n';
    }
    assert( isclose(sus_check[ii].mAngVel,sus[ii].mAngVel) );
    assert( isclose(sus_check[ii].mTorque,sus[ii].mTorque) );
  }
}

void testSusInteractions()
{
  std::cout << "Testing Poly Bacteria" << '\n';
  std::vector< RodShapedBacterium > sus;
  // std::vector< GenericInfected > inf;
  // std::vector< Phage >    phg;
  setUp(sus);
  testInteractions(sus);
}

// void testChainSusInteractions()
// {
//   std::cout << "Testing Chaining" << '\n';
//   std::vector< GenericInfected > inf;
//   std::vector< Phage >    phg;
//   std::vector< ChainingRodShapedBacterium > chaining;
//   setUp(chaining,inf,phg);
//   testInteractions(chaining);
// }

void testGeneralInteractions()
{
  SphericalBacterium sphere{0,0.8,0};
  RodShapedBacterium rod{0,0,0};

  std::cout << "rod: " << rod << '\n';
  // std::cout << "sphere: " << sphere << '\n';

  std::vector<IBacterium*> poly;
  poly.push_back(&sphere);
  poly.push_back(&rod);

  GridCells grid{ poly, constants::box_width };
  std::cout << "mGridSize "     << grid.mGridSize     << "\n"
            << "mNumGridCells " << grid.mNumGridCells << "\n"
            << "mCellSize "     << grid.mCellSize     << "\n"
            << "mOrigin "       << grid.mOrigin       << "\n";
  grid.updateGridCellLists(poly);
  grid.createNeighbourLists(poly);
  for ( auto &ii : poly )
  {
    for ( auto &nn : ii->getNeighbourList() )
    {
      std::cout << " ii: " << ii->getID()
                << " nn: " << nn->getID()
                << '\n';
    }
  }
  polyInteractParticles(
    poly
  );
  std::cout << "sphere: " << sphere.mForce << '\n';
  std::cout << "rod: " << rod.mForce << '\n';
  assert(isclose(sphere.mForce,-rod.mForce));
}

void testRodInteractions() {
  // Test a known interaction between two particles
  std::cout << "test pairCollision" << '\n';
  RodShapedBacterium A{0.0,-0.5,0.0,0.0,0.5*constants::pi,0,1.0};
  RodShapedBacterium B{0.0,1.5,0.0,0.0,0.5*constants::pi,0,2.0};
  // getForceAndTorqueOnCellA(&A,&B);
  // std::cout << "test1" << '\n';
  // std::cout << "force: " << A.mForce << '\n';
  // assert( isclose(max3(A.mForce,0.0),0.0) );
  // std::cout << "torque: " << A.mTorque << '\n';
  // assert( isclose(max3(A.mTorque,0.0),0.0) );
  //
  B.mAngles.x = 0.5*constants::pi;
  // std::cout << "test2" << '\n';
  // getForceAndTorqueOnCellA(&A,&B);
  // std::cout << "force: " << A.mForce << '\n';
  // assert( isclose(max3(A.mForce,0.0),0.0) );
  // std::cout << "torque: " << A.mTorque << '\n';
  // assert( isclose(max3(A.mTorque,0.0),0.0) );

  std::cout << "test3" << '\n';
  B.mLength=2.2;
  std::cout << A << '\n';
  std::cout << B << '\n';
  getForceAndTorqueOnCellA(&A,&B);

  double expected_force { 0.1*sqrt(0.1)/4 };
  std::cout << "force: " << A.mForce << " expected: " << expected_force << '\n';
  std::cout << "torque: " << A.mTorque << '\n';

  assert( isclose(0.0,fabs(A.mForce.x)) ); assert( isclose(0.0,fabs(A.mForce.z)) );
  assert( isclose(fabs(A.mForce.y),expected_force) );
  assert( isclose(max3(A.mTorque,0.0),0.0) );

  A.reset(); B.reset();
  std::cout << "test5" << '\n';
  B.mLength=2;
  A.mPos.x=-0.5*A.mLength-0.9; A.mPos.y=0.5;
  getForceAndTorqueOnCellA(&A,&B);
  std::cout << A.mForce << '\n';
  assert( isclose(0.0,fabs(A.mForce.y)) ); assert( isclose(0.0,fabs(A.mForce.z)) );
  assert( isclose(fabs(A.mForce.x),expected_force) );
  assert( isclose(max3(A.mTorque,0.0),0.0) );

  A.reset(); B.reset();
  getForceAndTorqueOnCellA(&B,&A);
  assert( isclose(0.0,fabs(B.mForce.y)) ); assert( isclose(0.0,fabs(B.mForce.z)) );
  std::cout << B.mForce << '\n';
  assert( isclose(fabs(B.mForce.x),expected_force) );
  assert( isclose(0.0,fabs(B.mTorque.y)) ); assert( isclose(0.0,fabs(B.mTorque.x)) );
  double expected_torque { expected_force };
  std::cout << "torque.z " << B.mTorque.z
            << " 0.5*length_B*force.x "
            << expected_torque
            << '\n';
  assert( isclose(B.mTorque.z , expected_torque ) );
}

int main(int argc, char const *argv[])
{
  // testRodInteractions();
  // testSusInteractions();
  testGeneralInteractions();
  return 0;
}
