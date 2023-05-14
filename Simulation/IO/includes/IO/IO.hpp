#ifndef IO_HPP
#define IO_HPP

// Standard libraries
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

// Custom
#include "RodShapedBacteria.hpp"
#include "MathUtility.hpp"
#include "constants.hpp"
#include "outputPaths.hpp"

extern std::vector<std::string> col_headers;
extern std::string sim_out_dir;

inline std::string createFileName(
  std::string stub,
  std::string stem="biofilm_",
  std::string dir =sim_out_dir
)
{
  return {
    dir + stem + stub + ".dat"
  };
}

inline std::string createFileName(
  uint output_count,
  std::string stem="biofilm_",
  std::string dir =sim_out_dir
)
{
  std::stringstream ss;
  ss << std::setw(5) << std::setfill('0') << output_count;
  return createFileName(ss.str(),stem,dir);
}

//--------------------------------- Outputs ------------------------------------
// template < class S >
// void printCellsToFile(
//   const std::string file_name,
//   const std::vector< S > &element_list,
//   bool verbose=false
// )
// {
//   // std::cout << "Saving to " << file_name << '\n';
//   std::ofstream out_file{ file_name, std::ios::out };
//   if (!out_file)
//   {
//     std::cerr << "Failed to open "
//               << file_name
//               << " for writing!"
//               << std::endl;
//     exit (10);
//   }
//   if ( verbose )
//     std::cout << "writing data to " << file_name << '\n';
//
//   for ( auto header : col_headers ) {
//     out_file << header;
//   }
//   out_file << "\n";
//
//   for ( auto &cell : element_list )
//   {
//     // std::cout << *cell << '\n';
//     printCellDataToFile(out_file,cell) << "\n";
//   }
//   out_file.close();
// }

template < class S >
void printCellsToFile(
  const std::string file_name,
  const std::vector< S* > &element_list,
  bool append=false,
  bool verbose=false
)
{
  // std::cout << "Saving to " << file_name << '\n';
  std::ofstream out_file;
  if ( append ) out_file.open( file_name, std::ios::app );
  else out_file.open( file_name, std::ios::out );

  if (!out_file)
  {
    std::cerr << "Failed to open "
              << file_name
              << " for writing!"
              << std::endl;
    exit (10);
  }
  if ( verbose )
    std::cout << "writing data to " << file_name << '\n';

  if ( !append )
    for ( auto header : col_headers ) {
      out_file << header;
    }
  // out_file << "\n";

  for ( auto &cell : element_list )
  {
    cell->printToFile(out_file);
  }
  out_file.close();
}

//---------------------------------- Inputs ------------------------------------

// template < class S >
void populateCellsFromFile(
  const std::string file_name,
  std::vector< IBacterium* > &element_list,
  bool verbose=false
);

/*--------------------------- Expected Types ---------------------------------*/

// template void populateCellsFromFile(
//   const std::string file_name,
//   std::vector<std::unique_ptr<RodShapedBacterium>> &element_list
// );
//
// template void printCellsToFile(
//   const std::string file_name,
//   const std::vector<std::unique_ptr<RodShapedBacterium>> &element_list
// );
//
// template void printCellsToFile(
//   const std::string file_name,
//   const std::vector<std::unique_ptr<GenericInfected>> &element_list
// );
//

/**

  Parameters:
    out: ostream to append to
    cell: pointer to cell to print out
  Effect:
    append cell data in the order:
      id, length, diameter, rcm, orientation
*/
template <class S>
std::ostream& printBaseCellDataToFile(std::ostream &out,
                                      const S &cell)
{
  // out << cell.mId               << "\t";
  out << cell.mLength           << "\t";
  out << cell.mRadius           << "\t";
  out << cell.mPos              << "\t";
  out << cell.getOrientation();
  return out;
}

template std::ostream& printBaseCellDataToFile(std::ostream &out,
                                               const RodShapedBacterium& cell);
// template std::ostream& printBaseCellDataToFile(std::ostream &out,
//                                                const ChainingRodShapedBacterium& cell);
// template std::ostream& printBaseCellDataToFile(std::ostream &out,
//                                                const GenericInfected& cell);
#endif
