// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <cassert>
#include <random>
#include <chrono>
#include <memory>

// Third party
#include "timer.hpp"

// Custom classes
#include "constants.hpp"         // definition of constants namespace
#include "MathUtility.hpp"
#include "RandUtil.hpp"
#include "RodShapedBacteria.hpp"
#include "SphericalBacteria.hpp"
// #include "VerletGrid.hpp"
// #include "biofilm_class.hpp"
// #include "forces.hpp"
// #include "test.hpp"

using vec = std::array<double,3>;
using Rvec = Vec3;

void testRodShapedBacteriumClass()
{

  // std::cout << "Checking basic utilities for RodShapedBacterium class" << '\n';
  RodShapedBacterium cellA{};
  // std::cout << "CoM at: " << cellA.mPos << '\n';
  // std::cout << "Direction at: " << cellA.getOrientation() << '\n';
  // std::cout << "cell area: " << cellA.getCellArea() << '\n';
  std::cout << "divide?: " << cellA.signalDivide() << " L " << cellA.mLength << '\n';
  cellA.grow(1.51e4);
  // std::cout << "cell area: " << cellA.getCellArea() << '\n';
  std::cout << "divide?: " << cellA.signalDivide() << " L " << cellA.mLength << '\n';

  // test basic movement
  cellA.mForce=Rvec{1,1,0};
  cellA.mTorque=Rvec{0,0,0.5};
  cellA.mAngVel=Vec3(0,0.5*constants::pi,0.0);
  cellA.mVel=Rvec{1,1,0};
  cellA.move(1);

  // test reset
  Rvec rcm { cellA.mPos };
  Rvec dir { cellA.getOrientation() };

  cellA.reset();
  cellA.move(0.1);
  Rvec rcm2 { cellA.mPos };
  Rvec dir2 { cellA.getOrientation() };
  assert ( static_cast<bool>( rcm2==rcm ) );
  assert ( static_cast<bool>( dir2==dir ) );

  std::cout << cellA.getPos() << '\n';

}
//
void testSphericalBacteriumClass()
{

  // std::cout << "Checking basic utilities for SphericalBacterium class" << '\n';
  SphericalBacterium cellA{2,27,5,0.1};
  IBacterium& iS = cellA;
  std::cout << "SphericalBacteria CoM at (getter): "
            << iS.getPos()
            << " attribute " << cellA.mPos << '\n';
  assert(cellA.mPos==iS.getPos());
  // std::cout << "Direction at: " << cellA.getOrientation() << '\n';
  // std::cout << "cell area: " << cellA.getCellArea() << '\n';
  std::cout << "divide?: " << cellA.signalDivide() << '\n';
  cellA.grow(10);
  // std::cout << "cell area: " << cellA.getCellArea() << '\n';
  // std::cout << "divide?: " << cellA.signalDivide() << '\n';

  // test basic movement
  // cellA.mForce=Rvec{1,1,0};
  // cellA.mVel=Rvec{1,1,0};
  // cellA.move(1);

  // test reset
  // Rvec rcm { cellA.mPos };
  // Rvec dir { cellA.getOrientation() };

  // cellA.reset();
  // cellA.move(0.1);
  // Rvec rcm2 { cellA.mPos };
  // Rvec dir2 { cellA.getOrientation() };
  // assert ( static_cast<bool>( rcm2==rcm ) );
  // assert ( static_cast<bool>( dir2==dir ) );

  // std::cout << cellA << '\n';

}

void testPoly()
{
  SphericalBacterium* sphere = new SphericalBacterium{2,27,5,0.1};
  RodShapedBacterium* rod = new RodShapedBacterium{30,18,101};
  rod->mLength=5.01;
  sphere->mRadius=1.78*2;

  std::vector<IBacterium*> access_vec;
  access_vec.push_back(sphere);
  access_vec.push_back(rod);

  std::cout << "testing positions" << '\n';
  for ( auto ii : access_vec )
  {
    std::cout << "I am " << ii->getMyType() << " " << ii->getPos() << '\n';
  }
  std::cout << "Expected: " << '\n';
  std::cout << sphere->mPos << '\n';
  std::cout << rod->mPos << '\n';

  std::cout << "testing divisions" << '\n';
  // Set this is to reverse order
  for ( int ii=access_vec.size()-1; ii>=0; --ii )
  {
    std::cout << "ii " << ii << " cells: " << access_vec.size() << '\n';
    std::cout << access_vec[ii]->getID()     << " "
              << access_vec[ii]->getMyType() << " "
              << access_vec[ii]->getPos()
              << '\n';
    if ( access_vec[ii]->signalDivide() )
    {
      access_vec[ii]->divide(access_vec);
    }
    std::cout << "access_vec: " << access_vec.size() << '\n';
  }
  for ( auto &cell : access_vec )
  {
    std::cout << "here" << '\n';
    std::cout << "cells: " << access_vec.size() << '\n';
    std::cout << cell->getID()     << " "
              << cell->getMyType() << " "
              << cell->getPos()
              << '\n';
    delete cell;
  }
}

int main(int argc, char const *argv[])
{
  RodShapedBacterium rodcell {};
  SphericalBacterium cell {};
  testSphericalBacteriumClass();
  testRodShapedBacteriumClass();
  testPoly();
  return 0;
}
