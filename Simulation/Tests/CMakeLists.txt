cmake_minimum_required (VERSION 3.16)

project (Tests LANGUAGES CXX VERSION 1.0)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(
 TestBac
 source/test_bacteria.cpp
 )
target_link_libraries(
 TestBac
 PRIVATE
 MathUtility
 Bacteria
)

add_executable(
 TestPolyGrid
 source/test_grid_poly.cpp
 )
target_link_libraries(
 TestPolyGrid
 PRIVATE
 MathUtility
 Grid
 Bacteria
)

add_executable(
 TestPolyInteractions
 source/testPolyInteractions.cpp
 )
target_link_libraries(
 TestPolyInteractions
 PRIVATE
 MathUtility
 Forces
 Grid
 Bacteria
 IO
)

add_executable(
  TestPolyBiofilm
  source/testPolyBiofilm.cpp
)
target_link_libraries(
  TestPolyBiofilm
  PRIVATE
  MathUtility
  Forces
  Bacteria
  Biofilm
  IO
)