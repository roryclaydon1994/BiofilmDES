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
// //#include <execution>

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

void initialiseBiofilm()
{
  // Check if there are any existing files in this directory, and if so begin
  // from here
  std::cout << "The sim directory is: " << sim_out_dir << '\n';

  int count { 0 };
  std::filesystem::path p1 { sim_out_dir };
  for (auto& p : std::filesystem::directory_iterator(p1)) ++count;

  std::stringstream ss;
  ss << std::setw(5) << std::setfill('0') << (count-2);

  std::string final_file { sim_out_dir+"final_"  +ss.str()+".dat" };
  std::string load_file  { sim_out_dir+"biofilm_"+ss.str()+".dat" };

  // Set the initial conditions
  std::vector<IBacterium*> access_vec;

  if ( std::filesystem::exists( final_file ) )
  {
    std::cout << "This directory " << final_file
              << " has already completed, exiting..." << '\n';
    exit(12);
  }
  else if ( std::filesystem::exists( load_file ) )
  {
    std::cout << "Attempting to reload from " << load_file << '\n';
    populateCellsFromFile(load_file,access_vec);
    exit(24);
  }
  else
  {
    std::cout << "No files found, restarting from a single cell." << '\n';
    access_vec.push_back(new RodShapedBacterium{0,0,0});
  }

  // SphericalBacterium* sphere = new SphericalBacterium{0,2,0};
  // access_vec.push_back(sphere);
  // access_vec.push_back(new RodShapedBacterium{2,0,0,});

#ifdef ADHESION
  // Adhesion via depletion is a longer range force and hence the box size is
  // increased accordingly
  const double interaction_length { 2*Ri+constants::nondim_avg_div_L };
  PolyBiofilm pb {
    access_vec,
    interaction_length
  };
#else // No adhesion
  PolyBiofilm pb {
    access_vec
  };
#endif // End ADHESION

// Currently not using this implementation but needs to be easier to integrate
// with the other simulation methods
#ifdef PHAGE
  // std::cout << "L: " << periodicBC(-0.5*constants::max_grid_size,constants::max_grid_size) << '\n';
  // const int num_phg { 40 };
  // std::vector<Phage*> phg_vec(num_phg,nullptr);
  // for (int ii=0;ii<num_phg;++ii)
  //   phg_vec[ii]=new Phage(-40,-20+ii,0);
  //
  // PolyBiofilm pb { access_vec,phg_vec };
#else // else no phage
#endif

  pb.runSim();
}

int main(int argc, char const *argv[])
{

#ifdef RANDOM_SEED
  std::cout << "setting random seed" << '\n';
  gen_rand.setSeed(
    std::chrono::high_resolution_clock::now().time_since_epoch().count()
  );
#endif

#if defined(ADHESION)
  std::cout << "Prevent running anything other than chaining" << '\n';
  exit(EXIT_FAILURE);
  if ( argc==6 )
  {
    const std::string run_dir    {            argv[1]  };
    const double kappa_depletion { std::stod( argv[2]) };
    const double Rd              { std::stod( argv[3]) };
    const double Ri              { std::stod( argv[4]) };
    const double rep_strength    { std::stod( argv[5]) };
    initialiseAdhesionParameters(kappa_depletion,Rd,Ri,rep_strength);
    sim_out_dir += "/" + run_dir + "/";
    initialiseBiofilm();
  }
  else
  {
    std::cout << "Expected 5 command line arguents! Received " << argc-1 << '\n';
    std::cout << "Example usage:\n"
              << "./main.out run_dir kappa_depletion Rd Ri rep_strength" << '\n';
    exit(EXIT_FAILURE);
  }
#elif defined(CHAINING)
  if ( argc==6 )
  {
    const std::string run_dir    {            argv[1]  }; // Run directory
    const double kappa           { std::stod( argv[2]) }; // Spring tension
    const double bend_rig        { std::stod( argv[3]) }; // Bending rigidity
    const double linking_prob    { std::stod( argv[4]) }; // Probability daughters link
    const double force_thresh    { std::stod( argv[5]) }; // Threshold force before breaking
    initialiseChainingParameters(kappa,bend_rig,linking_prob,force_thresh);
    sim_out_dir += "/" + run_dir + "/";
    initialiseBiofilm();
  }
  else
  {
    std::cout << "Expected 5 command line arguents! Received " << argc-1 << '\n';
    std::cout << "Example usage:\n"
              << "./main.out run_dir kappa bend_rig linking_prob force_thresh" << '\n';
    exit(EXIT_FAILURE);
  }
#elif defined(AG43)
  std::cout << "Prevent running anything other than chaining" << '\n';
  exit(EXIT_FAILURE);
  if ( argc==6 )
  {
    std::string run_dir    {            argv[1]  }; // Run directory
    double div_len         { std::stod( argv[2]) }; // division length
    double grwth_rate      { std::stod( argv[3]) }; // average growth rate
    double kappa           { std::stod( argv[4]) }; // Spring tension
    double force_thresh    { std::stod( argv[5]) }; // Threshold force before breaking
    initialiseAg43Parameters(div_len,grwth_rate,kappa,force_thresh);
    sim_out_dir += "/" + run_dir + "/";
    initialiseBiofilm();
  }
  else
  {
    std::cout << "Expected 5 command line arguents! Received " << argc-1 << '\n';
    std::cout << "Example usage:\n"
              << "./main.out run_dir divlen grwth_rate kappa force_thresh" << '\n';
    exit(1);
  }
#else
  std::cout << "Please compile a simulation mode: CHAINING, ADHESION, AG43" << '\n';
  exit(42);
#endif // End control input parameters

  return 0;
}
