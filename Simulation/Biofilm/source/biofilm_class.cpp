/*******************************************************************************
    BiofilmDES  - a program that simulates a growing colony of microbial cells

    Contributing author:
    Rory Claydon, University of Edinburgh, rory.claydon@ed.ac.uk

    Copyright (2020) The University of Edinburgh.

    The software is based on algorithms described in:

    Mechanically driven growth of quasi-two dimensional microbial colonies,
    F.D.C. Farrell, O. Hallatschek, D. Marenduzzo, B. Waclaw,
    Phys. Rev. Lett. 111, 168101 (2013).

    Three-dimensional distinct element simulation of spherocylinder crystallization.
    Pournin, L., Weber, M., Tsukahara, M. et al.
    Granul. Matter 7, 119–126 (2005).

    A fast algorithm to evaluate the shortest distance between rods,
    C. Vega, S. Lago,
    Comput. Chem., 18(1), 55-59 (1994)

    I would like to thank Bartlomiej Waclaw from Edinburgh University for some
    very useful discussions on algorithm stability, timestep choice and some
    potential optimisations to try out in future.

    This file is part of BiofilmDES.

    BiofilmDES is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BiofilmDES is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    A copy of the GNU General Public License can be found in the file
    License.txt or at <http://www.gnu.org/licenses/>.

    Compilation and run from current directory:
      make && ./biofilm.out 0 1.1 0.95

    Further details in the documentation

*******************************************************************************/

// Standard libraries
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <array>
#include <vector>
#include <utility>
#include <algorithm>
#include <cassert>
#include <random>
#include <memory>

// Third party
#include "timer.hpp"

// User defined
#include "VerletGrid.hpp"
#include "RodShapedBacteria.hpp"
#include "chaining_susceptible_class.hpp"
#include "biofilm_class.hpp"
#include "constants.hpp"        // definition of constants namespace
#include "MathUtility.hpp"        // definition of Vec3 class
#include "forces.hpp"     // definition of interactions
#include "IO.hpp"

template < class S >
Biofilm<S>::Biofilm (cellvec &proto_cells)
: mCellList{ proto_cells }
{
  // Expects only one type of bacterial cell for now, with all quantities in
  // it appropriately non-dimensionalised
  // mCellList.push_back(proto_cell);        // Add to the biofilm
}

// Delagate and only supply dt
template < class S >
Biofilm<S>::Biofilm ()
{
  mCellList.emplace_back();
}

template < class S >
Biofilm<S>::~Biofilm ()
{
  // Nothing required
}

template < class S >
double Biofilm<S>::size()
{
  double size { 0.0 };
  for ( auto& cell : mCellList )
  {
    size += cell.mLength + constants::nondim_rodSpheroDiam;
    // size += cell->mLength + cell->getDiameter();
  }
  return size;
}

template < class S >
void Biofilm<S>::updateOneTimeStep(double dt)
{

#ifdef DEBUG
  std::cout << "Grow" << '\n';
#endif
  std::for_each(
    std::execution::par_unseq,
    mCellList.begin(),
    mCellList.end(),
    [dt](S &cell){ cell.grow(dt); }
  );

#ifdef DEBUG
  std::cout << "CLL" << '\n';
#endif
  mTg.createCLL(mCellList);

#ifdef DEBUG
  std::cout << "Interactions" << '\n';
#endif
  interactParticles(
    mCellList,
    mTg.mCellBegin[IndexGrid::Types::SUS],
    mTg.mCellEnd[IndexGrid::Types::SUS]
  );

#ifdef DEBUG
  std::cout << "Move" << '\n';
#endif
  std::for_each(
    std::execution::par_unseq,
    mCellList.begin(),
    mCellList.end(),
    [dt](S &cell){ cell.move(dt); }
  );

#ifdef DEBUG
  std::cout << "Divide" << '\n';
#endif
  this->handleAllDivisionEvents();

}

// template < class S >
// long preallocCellList(Biofilm<S> &bio_film, double col_size_target)
// {
//   // auto& cell_list { bio_film.getCellsInFilm() };
//   // RodShapedBacterium& proto_cell { cell_list[0] };
//   // double avg_cell_len { 1 + 0.75*proto_cell.getDivLen() };
//   double avg_cell_len { 1 + 0.75*constants::nondim_avg_div_L };
//   double guess_elems { col_size_target / avg_cell_len } ;
//   long expected_element_number { static_cast<long>(ceil(guess_elems)) };
//   bio_film.getCellsInFilm().reserve(expected_element_number);
//   return expected_element_number;
// }

// This doesn't have to be a member function
template < class S >
void Biofilm<S>::evolveBiofilm(double col_size_target, long seed)
{

  // Seed the random number generator for the susceptible class
  Particle::mGenerator.seed(seed);

  // Find the expected number of elements and preallocate storage
  // long expected_element_number { preallocCellList(*this,col_size_target) };

  // Boxes must have at least a side length of the maximum distance between cells
  unsigned long file_index{0}; // output files sequentially labelled
  unsigned long ts{0};         // timesteps completed

  constexpr long output_timesteps {
    static_cast<long>( ceil( 0.02 / constants::dt ) )
  };
  std::cout << "OUT TS " << output_timesteps << '\n';
  const double dt { dC::dt*dC::rodModE/dC::zeta };
  std::cout << std::setprecision(10)
            << " sim dt: "  << dt
            << " real dt: " << dC::dt << " hrs" << '\n';

  while ( this->size() < col_size_target )
  {
    if ( ts % output_timesteps == 0 )
    {
      // Save data to file
      printCellsToFile(
        createFileName(file_index),
        mCellList,
        false
      );
      ++file_index;

      double col_size { this->size() };
      double progress { 100 * col_size / col_size_target };
      // TODO: Set a default
      std::cout << std::setprecision(3) << std::setw(5)
                << "\rcol_size: "    << col_size
                << " sim progress: " << progress
                << "% num elements: "<< mCellList.size()
                << " time (hrs): "   << ts*constants::dt << std::flush;
    }
    // Find all the forces and torques on all bacteria in this class
    // This updates the force and torque members of elements of mCellList
    try
    {
      this->updateOneTimeStep(dt);
    }
    catch (...)
    {
      // Save data to file
      printCellsToFile(
        createFileName("0","error_"),
        mCellList,
        true
      );
      exit(20);
    }
    // increment time steps
    ++ts;
  }
  // Save data to file
  printCellsToFile(
    createFileName("final"),
    mCellList,
    true
  );
}
//
// template < class S >
// void Biofilm<S>::printBiofilmStateToFile(const std::string file_name) const
// {
//   printCellsToFile<S>(file_name,mCellList);
// }
//
// template < class S >
// void printCellsToFile(const std::string file_name,
//                       const typename Biofilm<S>::cellvec &element_list)
// {
//   std::cout << "Saving to " << file_name << '\n';
//   std::ofstream biofilm_visfile{ file_name, std::ios::out };
//   if (!biofilm_visfile)
//   {
//       std::cerr << "Failed to open file for writing!" << std::endl;
//       exit (EXIT_FAILURE);
//   }
//
//   for ( auto header : S::col_headers ) {
//     biofilm_visfile << header;
//   }
//   biofilm_visfile << "\n";
//
//   for ( auto &cell : element_list )
//   {
//     //std::cout << *cell << '\n';
//     printCellDataToFile(biofilm_visfile,cell.get()) << "\n";
//   }
//   biofilm_visfile.close();
// }
//
// // Move this function to a static function for each type of susceptible and
// // use casting to reduce repetetion
// template <class RodShapedBacterium>
// void Biofilm<RodShapedBacterium>::createLogFile
// (
//   const std::string file_name,
//   const long seed
// )
// {
//   std::cout << "Saving logfile: " << file_name << '\n';
//   std::ofstream log_file{ file_name, std::ios::out };
//   if (!log_file)
//   {
//       std::cerr << "Failed to open file for writing!" << std::endl;
//       exit (EXIT_FAILURE);
//   }
//
//   log_file << "aspect\t";
//   log_file << "seed\n";
//   log_file << RodShapedBacterium::mAvgDivLen << "\t";
//   log_file << std::setprecision(25) << seed << '\n';
//
//   log_file.close();
//
// }
//
// // template <>
// // void Biofilm<ChainingRodShapedBacterium>::createLogFile
// // (
// //   const std::string file_name,
// //   const long seed
// // )
// // {
// //   std::cout << "Saving logfile: " << file_name << '\n';
// //   std::ofstream log_file{ file_name, std::ios::out };
// //   if (!log_file)
// //   {
// //       std::cerr << "Failed to open file for writing!" << std::endl;
// //       exit (EXIT_FAILURE);
// //   }
// //
// //   log_file << "aspect\t" << "kappa\t" << "bend_rig\t" << "link_prob\t" << "seed\n";
// //   log_file << RodShapedBacterium::mAvgDivLen      << '\t';
// //   log_file << ChainingRodShapedBacterium::mKappa       << '\t';
// //   log_file << ChainingRodShapedBacterium::mBendRig     << '\t';
// //   log_file << ChainingRodShapedBacterium::mLinkingProb << '\t';
// //   log_file << std::setprecision(25) << seed     << '\n';
// //   log_file.close();
// //
// // }
//
template < class S >
void Biofilm<S>::handleAllDivisionEvents()
{
  // Check cells which are about to divide (use reference to avoid copying large vector)
  // typename Biofilm<S>::cellvec& cells_in_film { bio_film.getCellsInFilm() };

   // store daughter cells to add to cells in film after checking for all divisions
  // daughter_cells.reserve(
  //   2 * ceil( 0.005 * static_cast<long>(cells_in_film.size()) )
  // );

  // Run backwards
  for ( int ii=static_cast<int>(mCellList.size()-1); ii>=0; --ii )
  {
    if ( mCellList[ii].signalDivide() )
    {
      // Handle iterator as the dividing cell will be removed from the list
      // This returns the iterator to the position after the removed element
      std::cout << "----------- " << ii << " -----------" << '\n';
      mCellList[ii].mId=ii;
      std::cout << mCellList[ii] << '\n';

      std::cout << "\n----------- Division result -----------" << '\n';
      createDaughters( ii, mCellList );
      std::cout << mCellList[0] << '\n';
      std::cout << mCellList[1] << '\n';

      // std::cout << *((it+1).base()) << '\n';
      // exit(1);
      // https://www.geeksforgeeks.org/how-to-erase-an-element-from-a-vector-using-erase-and-reverse_iterator/
      // mCellList.erase((it+1).base());
    }
  }
  // Add new daughter cells to the cell list
  // mCellList.insert
  // (
  //   mCellList.end(),
  //   daughter_cells.begin(),
  //   daughter_cells.end()
  // );
  // for ( auto &cell : daughter_cells )
  // {
  //   std::cout << cell.mAngles.x << " " << cell.mAngles.y << '\n';
  //   std::cout << cell.getOrientation() << '\n';
  //   assert( isclose(cell.getOrientation().z,0) );
  // }
}

// template < class S >
// typename void divideBacteria(
//   typename Biofilm<S>::cellvec &cell_list,
//   typename Biofilm<S>::cellvec &daughter_cells,
//   typename Biofilm<S>::cellvec::iterator it
// )
// {
//   // Make a reference to the mother cell
//   // RodShapedBacterium &mother_cell { *it };
//   createDaughters( (*it) , daughter_cells );
//
//   // Remove the mother from the environment
//   it = bio_film.getCellsInFilm().erase(it);
//
//   return it;
// }
