/*
  Test script

  Write functions to:
    - test diffusion - compare to MSD
    - test infection

  To compile:
    make

  To clean all compiled files:
    make clean

  To run:
    ./test.out
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

// User defined libraries
#include "constants.hpp"                    // definition of constants namespace
#include "MathUtility.hpp"                   // Vec3 defined here
#include "Phage.hpp"                  // Phage
#include "infected_class.hpp"               // GenericInfected
#include "RodShapedBacteria.hpp"               // RodShapedBacterium class
#include "IO_template.hpp"                  // Templates for I/O

void moveSimulationForward (Phage::phagevec &all_phage,
                           RodShapedBacterium::cellvec &all_susc,
                           GenericInfected::cellvec &all_inf,
                           double dt=0.1)
/*
Encapsulate moving forward in a function, such that simulating infection, lysis,
and time happens together.
*/
{
  //infect cells with phage inside (check default empty cellvec allinf works)
  handleAllInfectionEvents(all_phage,all_susc,all_inf);
  //Check whether any cells are ripe for lysing
  handleAllLysisEvents(all_phage, all_inf);

  //Introduce passage of time for infected
  for ( auto &infected : all_inf )
  {
    infected->mTimeSinceInfection += dt;
  }
}

void printSusToFile (const std::string file_name,
                     const RodShapedBacterium::cellvec &r_all_sus)

/*TODO Use template instead, for now used to print susceptibles to file using same
format as for phage and infected.*/

{
  std::ofstream sus_out{file_name, std::ios::out};
  if (!sus_out) {
      std::cerr << "Failed to open file for writing!" << std::endl;
      exit(EXIT_FAILURE);
  }
  // RC TODO: add to static variable
  for (int ii=0; ii < static_cast<int>(RodShapedBacterium::col_headers.size()); ++ii) {
    if (ii == static_cast<int>(RodShapedBacterium::col_headers.size()) - 1)
      sus_out << RodShapedBacterium::col_headers[ii];
    else
      sus_out << RodShapedBacterium::col_headers[ii] << "\t";
  }
  sus_out << "\n";

  for ( auto &susceptible : r_all_sus)
  {
    sus_out << *susceptible << "\n";
  }
  sus_out.close();
}

void testPhage(double dt = 0.1, long num_steps = 100, int num_phage = 10)
/*
test phage functionality, i.e multiple phage undergoing brownian motion*/
{

  Phage::phagevec all_phage;
  all_phage.reserve(num_phage); // preallocate memory

  // Make a bunch of randomly distributed phage
  std::uniform_real_distribution<double> rnd_pos(-10.0, 10.0);
  for ( int ii=0; ii<num_phage; ++ii )
  {
    all_phage.push_back(
      std::make_unique<Phage>(
        rnd_pos(Particle::mGenerator),
        rnd_pos(Particle::mGenerator),
        0
      )
    );
  }

  for ( long tt=0; tt < num_steps; ++tt)
  {

    //Save state of system every N timesteps
    if ( tt%10==0 )
    {
      std::stringstream ss;
      ss << std::setw(5) << std::setfill('0') << tt;
      std::string file_name{"./data/test_phage_motion_" + ss.str() +".dat"};
      printPhageToFile(file_name, all_phage);
    }

    // Find all the forces on all phage - mutates elements of phage list
    calcAllForces(all_phage,dt);

    // This could also be the same style as calcAllForces
    for ( auto &phage : all_phage )
    {
      phage->setVelFromForce(dt);
      phage->move(dt);
      phage->reset();
    }
  }

  // Make sure we print out the final timestep
  std::stringstream ss;
  ss << std::setw(5) << std::setfill('0') << num_steps;
  std::string file_name{"./data/test_phage_motion_" + ss.str() +".dat"};
  printPhageToFile(file_name, all_phage);
}

void testInfection(double dt = 0.1, long num_steps = 100)
/*
Phage infect bacteria
*/
{
  Phage::phagevec all_phage;
  RodShapedBacterium::cellvec all_susc;
  GenericInfected::cellvec all_inf;
  //Make 1 phage and bacteria that will infect
  all_phage.push_back(std::make_unique<Phage>(0,0,0));
  all_susc.push_back(std::make_unique<RodShapedBacterium>(0,0,0));

  //Make bacteria move forward in time
  for ( long tt=0; tt < num_steps; ++tt)
  {

    //Save state of system every N timesteps
    if ( tt%10==0 )
    {
      std::stringstream ss;
      ss << std::setw(5) << std::setfill('0') << tt;
      std::string file_name{"./data/test_infection_" + ss.str() +".dat"};
      printCellsToFile(file_name, all_inf);
    }
    //infect cells with phage inside (check default empty cellvec allinf works)
    handleAllInfectionEvents(all_phage,all_susc,all_inf);

    //Introduce passage of time for infected
    for ( auto &infected : all_inf )
    {
      infected->mTimeSinceInfection += dt;
    }
  }

}

void testLysis(double dt = 0.1, long num_steps = 100)
/*
infected bacteria undergo lysis and produce phage as expected.
*/
{
  Phage::phagevec all_phage;
  RodShapedBacterium::cellvec all_susc;
  GenericInfected::cellvec all_inf;
  //Make 1 phage and bacteria that will infect
  all_phage.push_back(std::make_unique<Phage>(0,0,0));
  all_susc.push_back(std::make_unique<RodShapedBacterium>(0,0,0));

  for ( long tt=0; tt < num_steps; ++tt)
  {

    //Save state of system every N timesteps
    if ( tt%10==0 )
    {

      std::stringstream ss;
      ss << std::setw(5) << std::setfill('0') << tt;

      //Print susceptible to test lysis
      std::string sus_file_name{"./data/test_bacteria_lysis_" + ss.str() +".dat"};
      printCellsToFile(sus_file_name, all_susc);

      //Print infected to test infection
      std::string inf_file_name{"./data/test_infection_" + ss.str() +".dat"};
      printCellsToFile(inf_file_name, all_inf);

      //Print phage to see if lytic cell produce them
      std::string ph_file_name{"./data/test_phage_production_" + ss.str() +".dat"};
      printPhageToFile(ph_file_name, all_phage);
    }

    //Move through a timestep (handle lys, inf, time)
    moveSimulationForward(all_phage, all_susc, all_inf, dt);

  }
}

void testIO()
{
  Phage::phagevec all_phage;
  RodShapedBacterium::cellvec all_susc;
  GenericInfected::cellvec all_inf;

  all_phage.push_back(std::make_unique<Phage>(2,3,4));
  all_susc.push_back(std::make_unique<RodShapedBacterium>(1,3,5));
  all_inf.push_back(std::make_unique<GenericInfected>(
    2,4,6,
    0,constants::pi*0.5,
    constants::nondim_rodGrwthRtePreFac,
    constants::nondim_init_length,
    100,4000));

  std::string sus_file_name{"./data/test_bacteria_IO.dat"};
  printCellsToFile(sus_file_name, all_susc);

  std::string inf_file_name{"./data/test_infected_IO.dat"};
  printCellsToFile(inf_file_name, all_inf);

  std::string ph_file_name{"./data/test_phage_IO.dat"};
  printPhageToFile(ph_file_name, all_phage);

}

int main(int argc, char const *argv[])
{
  testIO();
  double dt{0.1};
  long num_steps{50};
  testPhage(dt, num_steps, 10);
  testInfection(dt, num_steps);
  testLysis(dt, num_steps);
  return 0;
}
