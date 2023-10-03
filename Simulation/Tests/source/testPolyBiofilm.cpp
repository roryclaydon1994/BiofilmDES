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
#include "PolyBiofilm.hpp"

void testPolyBiofilm()
{

  // Set the initial conditions correctly for spherical bacteria
  SphericalBacterium* sphere = new SphericalBacterium{sqrt(2.0),0,0};
  RodShapedBacterium* rod = new RodShapedBacterium{0,0.5,0,constants::pi*0.25,constants::pi*0.5,1,2};

  std::vector<IBacterium*> initial_conditions;
  initial_conditions.push_back(rod);
  initial_conditions.push_back(sphere);

  PolyBiofilm pb { initial_conditions };
  pb.runSim();
}

int main(int argc, char const *argv[])
{
  testPolyBiofilm();
  return 0;
}
