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

#ifdef CHAINING

void createLogFile(
  std::string run_dir,
  uint ntrials,
  uint num_buckles
)
{
  std::filesystem::create_directories(run_dir);
  std::string filename=run_dir+'/'+"criticalLength.log";
  std::ofstream out_file;
  out_file.open( filename, std::ios::out );
  if (!out_file)
  {
    std::cerr << "Failed to open "
              << filename
              << " for writing!"
              << std::endl;
    exit (10);
  }
  std::vector<std::string> headers{};
  std::vector<double> values{};

  /*=== save basics ===*/
  headers.push_back("RodAspectRatio");
  values.push_back(constants::avg_div_L/constants::rodSpheroDiam);

  headers.push_back("RodModE");
  values.push_back(RodShapedBacterium::mRodModE);

  // In microns / hr
  headers.push_back("RodGrowthRate");
  values.push_back(constants::rodGrwthRtePreFac);

  headers.push_back("Kappa");
  values.push_back(RodShapedBacterium::mKappa);

  headers.push_back("BendRig");
  values.push_back(RodShapedBacterium::mBendRig);

  headers.push_back("LinkingProb");
  values.push_back(RodShapedBacterium::mLinkingProb);

  headers.push_back("ForceThresh");
  values.push_back(RodShapedBacterium::mForceThresh);

  headers.push_back("NTrials");
  values.push_back(ntrials);

  headers.push_back("NBuckles");
  values.push_back(num_buckles);

  for ( uint ii=0; ii<headers.size(); ++ii )
  {
    if ( ii<headers.size()-1 )
      out_file << headers[ii] << '\t';
    else
      out_file << headers[ii] << '\n';
  }

  for ( uint ii=0; ii<values.size(); ++ii )
  {
    if ( ii<headers.size()-1 )
      out_file << values[ii] << '\t';
    else
      out_file << values[ii] << '\n';
  }
  out_file.close();
  std::cout << "\nSaved log file to " << filename << '\n';
}


void saveLengths(
  std::string run_dir,
  const std::vector<double> &lengths
)
{
  std::string filename=run_dir+'/'+"lengths.dat";
  std::ofstream out_file;
  out_file.open( filename, std::ios::out );
  if (!out_file)
  {
    std::cerr << "Failed to open "
              << filename
              << " for writing!"
              << std::endl;
    exit (10);
  }
  for ( uint ii=0; ii<lengths.size(); ++ii )
  {
    if ( ii<lengths.size()-1 )
      out_file << lengths[ii] << '\t';
    else
      out_file << lengths[ii] << '\n';
  }
  out_file.close();
  std::cout << "Saved log file to " << filename << '\n';
}

void saveTimes(
  std::string run_dir,
  const std::vector<double> &times
)
{
  std::string filename=run_dir+'/'+"times.dat";
  std::ofstream out_file;
  out_file.open( filename, std::ios::out );
  if (!out_file)
  {
    std::cerr << "Failed to open "
              << filename
              << " for writing!"
              << std::endl;
    exit (10);
  }
  out_file << "Critical times (hrs)" << '\n';
  for ( uint ii=0; ii<times.size(); ++ii )
  {
    if ( ii<times.size()-1 )
      out_file << times[ii] << '\t';
    else
      out_file << times[ii] << '\n';
  }
  out_file.close();
  std::cout << "Saved log file to " << filename << '\n';
}

uint getNumChains(const std::vector<IBacterium*> &cells)
{
  uint num_chains { 0 };
  for ( auto &cell : cells )
  {
    // Count the number of chain heads
    if ( cell->getUpperLink()!=nullptr && cell->getLowerLink()==nullptr )
      ++num_chains;
  }
  return num_chains;
}

uint getNumSingles(const std::vector<IBacterium*> &cells)
{
  uint num_singles { 0 };
  for ( auto &cell : cells )
  {
    // Count the number of chain heads
    if ( cell->getUpperLink()==nullptr && cell->getLowerLink()==nullptr )
      ++num_singles;
  }
  return num_singles;
}

bool checkBuckled(const std::vector<IBacterium*> &cells)
{
  bool has_buckled { false };
  for ( uint ii=0; ii<cells.size()-1; ++ii )
  {
    const Vec3 t1 { cells[ii  ]->getOrientation() };
    const Vec3 t2 { cells[ii+1]->getOrientation() };
    // Previously was 5 degrees
    if ( acos(dot(t1,t2))*180/constants::pi > 3 )
    {
      has_buckled=true;
      break;
    }
  }
  return has_buckled;
}

void runBucklingSim(
  uint &num_buckles,
  std::vector<double> &critical_lengths,
  std::vector<double> &critical_times,
  const uint trial_num,
  const std::string save_dir
)
{
  // Set the initial conditions
  std::vector<IBacterium*> access_vec;
  RodShapedBacterium* rod = new RodShapedBacterium{0,0,0};
  access_vec.push_back(rod);

  // Set up biofilm
  PolyBiofilm pb {
    access_vec
  };

  uint verlet_counter { 0 };      // Count successful steps since last rebinning
  bool update_neighbours { true };// Always bin neighbours on first step

  bool has_buckled { false };
  uint num_chains { 0 };

  uint counter { 0 };
  uint out_counter { 0 };
  while ( has_buckled==false && num_chains<2 )
  {
    if ( counter%pb.mOutFreq==0 )
    {
      std::cout << "\rsaving at time: "
                << std::left << std::setw(4) << std::setprecision(3)
                << counter*constants::dt << "(h) "
                << "num bacteria: " << std::setw(6)
                << std::left << pb.mCells.size()
                << std::left << " num chains: " << num_chains
                << " outs: " << out_counter
                << std::flush;
      ++out_counter;
    }

    // refresh list if every N timesteps
    if ( verlet_counter>=pb.mGrid.mVerletUptime )
    {
      // std::cout << "Binning: verlet_counter>=mGrid.mVerletUptime" << '\n';
      update_neighbours=true;
    }
    else ++verlet_counter;

    pb.updateOneTimeStep(update_neighbours,verlet_counter);

    num_chains = getNumChains(pb.mCells);
    if ( num_chains==1 ) has_buckled = checkBuckled(pb.mCells);
    if ( pb.mCells.size()>1 && getNumSingles(pb.mCells)>0 ) break;

    ++counter;
  }
  if ( has_buckled )
  {
    ++num_buckles;
    double chain_length { 0.0 };
    for ( auto cell : pb.mCells )
      chain_length+=cell->getLength()+2*cell->getRadius();
    critical_lengths.push_back(chain_length);
    critical_times.push_back(counter*constants::dt);
  }

  std::stringstream ss;
  ss << std::boolalpha << has_buckled
     << std::setw(5) << std::setfill('0') << trial_num;
  printCellsToFile(
    createFileName(
      out_counter,
      "buckled_"+ss.str()+"_",
      save_dir
    ),
    pb.mCells,
    false,
    false);
}

void runCriticalLengthTrials(
  const double kappa,
  const double bend_rig,
  const double force_thresh
)
{
  constexpr double linking_sep { 0.5e-2 };
  double linking_prob { 1.0 };

  initialiseChainingParameters(
    kappa,
    bend_rig,
    linking_prob,
    force_thresh
  );

  constexpr uint n_trials { 200 };
  uint counter { 0 };
  while ( RodShapedBacterium::mLinkingProb>=0.0 )
  {
    std::cout << "Running for lp: " << RodShapedBacterium::mLinkingProb << '\n';

    std::stringstream ss;
    ss << std::setw(5) << std::setfill('0') << counter;
    std::string run_dir { sim_out_dir + "/trial_" + ss.str() + '/' };
    std::filesystem::create_directories(run_dir);

    uint num_buckles { 0 };
    std::vector<double> critical_lengths;
    critical_lengths.reserve(n_trials);

    std::vector<double> critical_times;
    critical_times.reserve(n_trials);
    for ( uint ii=0; ii<n_trials; ++ii )
    {
      std::cout << "trial: " << ii << '\n';
      gen_rand.setRandomSeed();
      runBucklingSim(num_buckles,critical_lengths,critical_times,ii,run_dir);
    }

    createLogFile(run_dir,n_trials,num_buckles);
    saveLengths(run_dir,critical_lengths);
    saveTimes(run_dir,critical_times);

    // // If this chain has no chance of buckling do not check lower lps
    if ( num_buckles == 0 )
      break;

    RodShapedBacterium::mLinkingProb-=linking_sep;
    ++counter;
  }
}

#endif

int main(int argc, char const *argv[])
{
#ifndef CHAINING
  std::cout << "Please recompile " << '\n';
  exit(2);
#else
  if ( argc==5 )
  {
    const std::string run_dir    {            argv[1]  }; // Run directory
    const double kappa           { std::stod( argv[2]) }; // Spring tension
    const double bend_rig        { std::stod( argv[3]) }; // Bending rigidity
    const double force_thresh    { std::stod( argv[4]) }; // Threshold force before breaking
    sim_out_dir += "/" + run_dir + "/";
    runCriticalLengthTrials(kappa,bend_rig,force_thresh);
  }
  else
  {
    std::cout << "Expected 4 command line arguents! Received " << argc-1 << '\n';
    std::cout << "Example usage:\n"
              << "./main.out run_dir kappa bend_rig force_thresh" << '\n';
    exit(EXIT_FAILURE);
  }

  return 0;
#endif
}
