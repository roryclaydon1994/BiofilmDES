cmake_minimum_required (VERSION 3.16)

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR
    "You cannot build in a source directory (or any directory with a CMakeLists.txt file).
    Please make a build subdirectory.
    Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# set the project name and version
if (CMAKE_VERSION VERSION_LESS 3.0)
  #cmake_policy(SET CMP0048 OLD)
  project (Tests )
else()
  cmake_policy(SET CMP0048 NEW)
  project (Tests LANGUAGES CXX VERSION 1.0)
endif()


# specify the C++ standard
#set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD_REQUIRED True)

#add_executable(
#  TestBac
#  source/test_bacteria.cpp
#  )
#target_link_libraries(
#  TestBac
#  PRIVATE
#  MathUtility
#  Bacteria
#)

#add_executable(
#  TestPolyGrid
#  source/test_grid_poly.cpp
#  )
#target_link_libraries(
#  TestPolyGrid
#  PRIVATE
#  MathUtility
#  Grid
#  Bacteria
#)

#add_executable(
#  TestGrid
#  source/test_grid_poly.cpp
#  )
#target_link_libraries(
#  TestGrid
#  PRIVATE
#  MathUtility
#  Grid
#  Bacteria
#  Phage
#)
#
#add_executable(
#  TestInteractions
#  source/test_forces.cpp
#  )
#target_link_libraries(
#  TestInteractions
#  PRIVATE
#  MathUtility
#  Grid
#  Bacteria
#  IO
#  Phage
#  tbb
#)
#


#add_executable(
#  TestPolyInteractions
#  source/testPolyInteractions.cpp
#  )
#target_link_libraries(
#  TestPolyInteractions
#  PRIVATE
#  MathUtility
#  Grid
#  Forces
#  Bacteria
#  IO
#  Phage
#  tbb
#)

add_executable(
  TestPolyBiofilm
  source/testPolyBiofilm.cpp
)
target_link_libraries(
  TestPolyBiofilm
  PRIVATE
  MathUtility
  Forces
  Phage
  Bacteria
  Biofilm
  IO
)

#
#add_executable(
#  TestBiofilm
#  source/testBiofilm.cpp
#)
#target_link_libraries(
#  TestBiofilm
#  PRIVATE
#  MathUtility
#  Particle
#  Phage
#  Grid
#  Bacteria
#  Biofilm
#  IO
#)

#add_executable(
#  TestInf
#  source/test_infection.cpp
#  )
#target_link_libraries(TestInf PRIVATE MathUtility Particle Phage Grid Bacteria Biofilm ${Boost_LIBRARIES})
