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
#include <omp.h>
#include <filesystem>

// Custom classes
#include "constants.hpp"         // definition of constants namespace
#include "MathUtility.hpp"
#include "RandUtil.hpp"
#include "IBacterium.hpp"
#include "SphericalBacteria.hpp"
#include "RodShapedBacteria.hpp"
#include "VerletGrid.hpp"
#include "forces.hpp"
#include "IO.hpp"
#include "PolyBiofilm.hpp"

std::vector<IBacterium*> initialiseBiofilm()
{
  // Set the initial conditions
  std::vector<IBacterium*> initial_condtions;

  std::cout << "Starting from a single cell." << '\n';
  auto* rod = new RodShapedBacterium{
    0,0,0,
    0,constants::pi*0.5,
    RodShapedBacterium::mAvgGrwthRate,
    0.5*(RodShapedBacterium::mAvgDivLen-2*RodShapedBacterium::mRadius)
  };
  initial_condtions.push_back(rod);

  return initial_condtions;
}

int main(int argc, char const *argv[])
{

#ifdef RANDOM_SEED
  std::cout << "setting random seed" << '\n';
  gen_rand.setSeed(
    std::chrono::high_resolution_clock::now().time_since_epoch().count()
  );
#endif

#if defined(CHAINING)
  std::string run_dir; // Run directory
  if ( argc==6 )
  {
    run_dir = argv[1];                              // Run directory
    double kappa           { std::stod( argv[2]) }; // Spring tension
    double bend_rig        { std::stod( argv[3]) }; // Bending rigidity
    double linking_prob    { std::stod( argv[4]) }; // Probability daughters link
    double force_thresh    { std::stod( argv[5]) }; // Threshold force before breaking
    initialiseChainingParameters(kappa,bend_rig,linking_prob,force_thresh);
  }
  else if ( argc==8 )
  {
    run_dir = argv[1];                              // Run directory
    double kappa           { std::stod( argv[2]) }; // Spring tension
    double bend_rig        { std::stod( argv[3]) }; // Bending rigidity
    double linking_prob    { std::stod( argv[4]) }; // Probability daughters link
    double force_thresh    { std::stod( argv[5]) }; // Threshold force before breaking
    double aspect_ratio    { std::stod( argv[6]) }; // Division length of the bacterium
    double growth_rate     { std::stod( argv[7]) }; // Growth rate

    initialiseRodParameters(aspect_ratio,growth_rate);
    initialiseChainingParameters(kappa,bend_rig,linking_prob,force_thresh);
  }
  else
  {
    std::cout << "Expected 5 command line arguents! Received " << argc-1 << '\n';
    std::cout << "Example usage:\n"
              << "./main.out run_dir kappa bend_rig linking_prob force_thresh" << '\n';
    exit(EXIT_FAILURE);
  }
  if ( !std::filesystem::exists(sim_out_dir) )
  {
    std::filesystem::create_directories(sim_out_dir);
  }

  sim_out_dir += "/" + run_dir + "/";
  std::vector<IBacterium*> initial_conditions{ initialiseBiofilm() };
  PolyBiofilm pb { initial_conditions };
  pb.runSim();
#else
  std::cout << "Please define one of the following MACROS: CHAINING" << '\n';
  exit(42);
#endif // End control input parameters

  return 0;
}
