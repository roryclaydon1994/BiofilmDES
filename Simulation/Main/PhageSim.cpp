/*
  Description
  ./main.out dt total_time lysis_period burst_size phage_diffusion
  ./main.out 1e-4 10 5 5 0.02
*/

// Standard libraries
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <memory>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <chrono>

// User defined libraries
#include "constants.hpp"                    // definition of constants namespace
#include "MathUtility.hpp"                   // Vec3 defined here
#include "Phage.hpp"                  // Phage
#include "infected_class.hpp"               // GenericInfected
#include "RodShapedBacteria.hpp"               // RodShapedBacterium class
#include "IO_template.hpp"                  // Templates for I/O
#include "infected_biofilm_class.hpp"       // InfectedBiofilm

// This should also be in a separate file
void runSim(std::string sus_file, double dt, long timesteps,
            long output_frequency, double lysis_period, int burst_size,
            double phage_diffusion, char centering,
            std::string base_output_path="")
{

  // Need to move this to biofilm_IO later
  InfectedBiofilm biofilm {};
  populateCellsFromFile(sus_file,biofilm.mSusList);

  if (centering == 's')
  {
    // Find the furthest deviation from the origin mapped to the first quadrant
    auto X_max = std::max_element(
      biofilm.mSusList.begin(),biofilm.mSusList.end(),
        []( const auto &particleA, const auto &particleB )
          {
            return particleA->mPos.mX < particleB->mPos.mX;
          } );
    Vec3 init_phg_rcm { (*X_max)->mPos };
    biofilm.mPhgList.push_back(std::make_unique<Phage>(
    init_phg_rcm, lysis_period, burst_size, phage_diffusion));
  }

  else
  {
    Vec3 init_phg_rcm { 0,0,0 };
    biofilm.mPhgList.push_back(std::make_unique<Phage>(
    init_phg_rcm, lysis_period, burst_size, phage_diffusion));
  }

  handleAllInfectionEvents(biofilm.mPhgList,biofilm.mSusList,biofilm.mInfList);

  std::cout << "------------------------------------------------------" << '\n';
  std::cout << "------------------Begin Simulation--------------------" << '\n';
  std::cout << "------------------------------------------------------" << '\n';

  if ( output_frequency > timesteps )
  {
    std::cout << "Error - output step size less than time_steps!" << '\n';
    std::cout << "Output frequency: " << output_frequency
              << " timesteps: " << timesteps
              << '\n';
    std::cout << "Exiting..." << '\n';
    exit (EXIT_FAILURE);
  }
  std::cout << "outputting every " << output_frequency*dt
            << " or every " << output_frequency*dt*(constants::zeta/constants::rodModE)
            << " hours "  << '\n';
  long output_counter { 0 };
  const long num_outputs { static_cast<long>(ceil(timesteps/output_frequency)) };
  for ( long tt=0; tt<=timesteps; ++tt )
  {
    // check for outputting
    if ( tt%output_frequency==0 )
    {
      std::cout << "\rsaving output: "
                << std::setw(5) << std::right << output_counter
                << " / "
                << std::left << std::setw(5) << num_outputs
                << " ---- "
                << "phg: " << std::setw(7) << std::left << biofilm.mPhgList.size() << " "
                << "inf: " << std::setw(7) << std::left << biofilm.mInfList.size() << " "
                << "sus: " << std::setw(7) << std::left << biofilm.mSusList.size()
                << std::flush;
      std::stringstream ss;
      ss << std::setw(5) << std::setfill('0') << output_counter;

      // RC: This directory should be passed as a command line argument
      // std::string base_output_path =
      // "./results/init_large_#.txt-0.0025-400000-10000-5-1-c_lowres/Colony0/";

      //Print susceptible to test lysis
      std::string sus_file_name {
        base_output_path + "main_bacteria_" + ss.str() +".dat"
      };
      printCellsToFile(sus_file_name, biofilm.mSusList);

      //Print infected to test infection
      std::string inf_file_name {
        base_output_path + "main_infected_" + ss.str() +".dat"
      };
      printCellsToFile(inf_file_name, biofilm.mInfList);

      //Print phage to see if lytic cell produce them
      std::string ph_file_name {
        base_output_path + "main_phage_" + ss.str() +".dat"
      };
      printPhageToFile(ph_file_name, biofilm.mPhgList);

      ++output_counter;

      if ( tt==timesteps )
      {
        break;
      }

      if (
        ( static_cast<long>(biofilm.mSusList.size()) == 0 )
        &&
        ( static_cast<long>(biofilm.mInfList.size()) == 0 )
      )
      {
        std::cout << "Biofilm wiped out!" << '\n';
        break;
      }

      if (
        ( static_cast<long>(biofilm.mInfList.size()) == 0 )
        &&
        ( static_cast<long>(biofilm.mPhgList.size()) == 0 )
      )
      {
        std::cout << "Biofilm survived!" << '\n';
        break;
      }
    }

    // update biofilm one dt
    biofilm.updateOneTimeStep(dt);

  }
  std::cout << '\n';
}

int main(int argc, char const *argv[]) {
  if ( argc != 9 )
  {
    std::cout << "Usage: " << argv[0] << " ";
    std::cout << "sus_file ";
    std::cout << "dt ";
    std::cout << "total_time ";
    std::cout << "lysis_period ";
    std::cout << "burst_size ";
    std::cout << "phage_diffusion ";
    std::cout << "centering"<< '\n';
    std::cout << "Example: "
              << argv[0]
              << " ./initial_conditions/init_large_0.txt "
              << "2.5e-3 250000 10000 5 1 c ./results/" << '\n';
    exit (EXIT_FAILURE);
  }
  std::string sus_file   { argv[1] };
  double dt              { std::stod(argv[2]) };
  double total_time      { std::stod(argv[3]) };

  // RC: TODO: shift these to static variables in the phage class later
  double lysis_period    { std::stod(argv[4]) };
  double burst_size      { std::stod(argv[5]) };
  double phage_diffusion { std::stod(argv[6]) };
  char centering         { *argv[7] };
  std::string out_dir    {  argv[8] };

  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  Particle::mGenerator.seed( seed );

  std::cout << "dt: " << dt << " is "
            << dt*(constants::zeta/constants::rodModE) << " hours" << '\n';
  std::cout << "total_time: " << total_time << " is "
            << total_time*constants::zeta/constants::rodModE << " hours" << '\n';
  std::cout << "lysis period: " << lysis_period << " is "
            << lysis_period*constants::zeta/constants::rodModE << " hours" << '\n';
  std::cout << "diff_const in: " << phage_diffusion << " is "
            <<  phage_diffusion * pow(constants::rodSpheroDiam,2)
              /  ( 3600 * ( constants::zeta/constants::rodModE ) )
            << " in microns^2/s " << "\n";

  long expected_lysis_generations { static_cast<long>(total_time/lysis_period) };
  if ( expected_lysis_generations <= 10 )
  {
    std::cout << "\n-----------------------Warning----------------------\n";
    std::cout << "-------- Only "
              << expected_lysis_generations
              << " expected_lysis_generations ---------" << '\n';
    std::cout << "consider increasing total_time to at least "
              << lysis_period*10 << '\n';
    std::cout << "----------------------------------------------------" << "\n";
  }

  // Characteristic time is smallest of either average time taken to lyse or time
  // taken for a phage to move some \Delta x
  // (to avoid missing many infections or flying through bacteria)
  double stdev_phage_step { 0.5 }; // Standard deviation of random walk in microns
  const double shortest_characteristic_time {
    std::min(lysis_period,pow(stdev_phage_step,2)/phage_diffusion)
  };

  if ( dt<1e-3*shortest_characteristic_time )
  {
    std::cout << "\n-----------------------Warning----------------------\n";
    std::cout << "dt " << dt
              << " is quite small, consider increasing to "
              << 1e-2*shortest_characteristic_time << '\n';
    std::cout << "----------------------------------------------------" << "\n";
  }
  if ( dt>1e-1*shortest_characteristic_time )
  {
    std::cout << "\n-----------------------Warning----------------------\n";
    std::cout << "dt " << dt
              << " is very large, consider decreasing to "
              << 1e-2*shortest_characteristic_time << '\n';
    std::cout << "----------------------------------------------------" << "\n";
  }

  const double output_every_hrs {
    // std::min(
    //   1e-2*constants::rodModE/constants::zeta,
    //   1e-1*shortest_characteristic_time
    // )
    // 1e-1*shortest_characteristic_time
    0.1*lysis_period
  };

  long output_frequency {
    static_cast<long>( ceil(output_every_hrs/dt) )
  };
  if ( output_frequency<=10 )
  {
    std::cout << "----------------------- Error ----------------------" << '\n';
    std::cout << "Check input params! Output frequency is "
              << output_frequency << '\n';
    std::cout << "Exiting..." << '\n';
    std::cout << "----------------------------------------------------" << "\n\n";
    exit(EXIT_FAILURE);
  }
  if ( centering != 'c' && centering != 's' )
  {
    std::cout << "----------------------- Error ----------------------" << '\n';
    std::cout << "Infection centering input not valid! Please choose from 'c' for"
              << " center or 's' for side";
    std::cout << "Exiting..." << '\n';
    std::cout << "----------------------------------------------------" << "\n\n";
    exit(EXIT_FAILURE);
  }

  long timesteps { static_cast<long>( ceil( total_time/dt ) ) };
  const long num_outputs { static_cast<long>(ceil(timesteps/output_frequency)) };
  if ( num_outputs >= 1e6 )
  {
    std::cout << "Too many outputs requested: " << num_outputs << '\n';
    exit(EXIT_FAILURE);
  }

  /*------------------------ Create log file ---------------------------------*/
  std::string filename { out_dir+"infected_biofilm.log" };
  std::ofstream out_file{ filename, std::ios::out };
  if (!out_file)
  {
      std::cerr << "Failed to open file for writing!" << std::endl;
      exit (EXIT_FAILURE);
  }
  out_file << "sus_file"        << '\t';
  out_file << "dt"              << '\t';
  out_file << "total_time"      << '\t';
  out_file << "lysis_period"    << '\t';
  out_file << "burst_size"      << '\t';
  out_file << "phage_diffusion" << '\t';
  out_file << "centering"       << '\t';
  out_file << "seed"            << '\n';
  out_file << sus_file          << '\t';
  out_file << dt                << '\t';
  out_file << total_time        << '\t';
  out_file << lysis_period      << '\t';
  out_file << burst_size        << '\t';
  out_file << phage_diffusion   << '\t';
  out_file << centering         << '\t';
  out_file << seed              << '\n';
  out_file.close();

  runSim(
    argv[1], dt, timesteps, output_frequency,
    lysis_period, burst_size, phage_diffusion, centering,
    out_dir
  );

  return 0;
}
