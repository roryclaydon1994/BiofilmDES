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
#include "IBacterium.hpp"
#include "SphericalBacteria.hpp"
#include "RodShapedBacteria.hpp"
#include "VerletGrid.hpp"
#include "forces.hpp"
#include "IO.hpp"
#include "PolyBiofilm.hpp"

// void setUp()
// {
//   // Read in from file
//   populateCellsFromFile(
//     "/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/InitialConditions/init_large_1.txt",
//     sus,
//     true
//   );
//   // Set the grid and update the HyperParams
//   tg.createCLL(sus,inf,phg);
//   return;
// }

void testPolyBiofilm()
{

  // Set the initial conditions correctly for spherical bacteria
  // access_vec.push_back(new RodShapedBacterium{0,0.5,0,constants::pi*0.25,constants::pi*0.5,1,2});
  SphericalBacterium* sphere = new SphericalBacterium{sqrt(2.0),0,0};
  RodShapedBacterium* rod = new RodShapedBacterium{0,0.5,0,constants::pi*0.25,constants::pi*0.5,1,2};

  std::vector<IBacterium*> access_vec;
  access_vec.push_back(rod);
  access_vec.push_back(sphere);
  // access_vec.push_back(new SphericalBacterium{0,0,sqrt(2.0)});
  // access_vec.push_back(new RodShapedBacterium{-2,0,0});

#ifdef PHAGE
  // std::cout << "L: " << periodicBC(-0.5*constants::max_grid_size,constants::max_grid_size) << '\n';
  // const int num_phg { 40 };
  // std::vector<Phage*> phg_vec(num_phg,nullptr);
  // for (int ii=0;ii<num_phg;++ii)
  //   phg_vec[ii]=new Phage(-40,-20+ii,0);
  //
  // PolyBiofilm pb { access_vec,phg_vec };
  PolyBiofilm pb { access_vec };

#else
  PolyBiofilm pb { access_vec };
#endif
  pb.runSim();
}

int main(int argc, char const *argv[])
{
  testPolyBiofilm();
  return 0;
}
