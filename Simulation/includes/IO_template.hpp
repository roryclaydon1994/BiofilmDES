#ifndef IO_TEMPLATE
#define IO_TEMPLATE

// Standard libraries
#include <memory>
#include <vector>

template < class S >
void printCellsToFile(
  const std::string file_name,
  const std::vector< std::unique_ptr< S >  > &element_list
);

template < class S >
void populateCellsFromFile(
  const std::string sus_file_name,
  std::vector< std::unique_ptr< S >  > &element_list
);
/**<
    Parameters:
    @param[in]  file_name: file for which to load data from
    @param[out] mSusList,mInfList,mPhgList to match cells in the file.
    @returns Void
*/

#endif
