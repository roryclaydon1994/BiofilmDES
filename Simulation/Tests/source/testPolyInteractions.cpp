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

// Custom classes
#include "constants.hpp"         // definition of constants namespace
#include "MathUtility.hpp"
#include "IBacterium.hpp"
#include "SphericalBacteria.hpp"
#include "RodShapedBacteria.hpp"
#include "VerletGrid.hpp"
#include "forces.hpp"
#include "IO.hpp"

void testGeneralInteractions()
{
  SphericalBacterium sphere{0,0.8,0};
  RodShapedBacterium rod{0,0,0};

  std::cout << "rod: " << rod << '\n';

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
  B.mAngles.x = 0.5*constants::pi;

  std::cout << "test3" << '\n';
  B.mLength=2.2;
  std::cout << A << '\n';
  std::cout << B << '\n';

  RodShapedBacterium* pA { &A };
  const RodShapedBacterium* pB { &B };
  pairCollision(pA,pB);

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
  pairCollision(&A,&B);
  std::cout << A.mForce << '\n';
  assert( isclose(0.0,fabs(A.mForce.y)) ); assert( isclose(0.0,fabs(A.mForce.z)) );
  assert( isclose(fabs(A.mForce.x),expected_force) );
  assert( isclose(max3(A.mTorque,0.0),0.0) );

  A.reset(); B.reset();
  pairCollision(&B,&A);
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
  testGeneralInteractions();
  testRodInteractions();
  return 0;
}
