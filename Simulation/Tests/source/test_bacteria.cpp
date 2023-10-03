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

using vec = std::array<double,3>;
using Rvec = Vec3;

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
  testPoly();
  return 0;
}
