cmake_minimum_required (VERSION 3.16)

project (Bacteria LANGUAGES CXX VERSION 1.0)

add_library(
  Bacteria STATIC
  source/RodShapedBacteria.cpp
  source/SphericalBacteria.cpp
  source/Divisions.cpp
  includes/Bacteria/RodShapedBacteria.hpp
  includes/Bacteria/SphericalBacteria.hpp
  )
target_include_directories(Bacteria PUBLIC includes/Bacteria)
target_link_libraries(Bacteria PRIVATE MathUtility IO Phage)
