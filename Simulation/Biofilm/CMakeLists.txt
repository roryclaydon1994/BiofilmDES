cmake_minimum_required (VERSION 3.16)

add_library(
  Biofilm
  STATIC
  includes/Biofilm/PolyBiofilm.hpp
  source/PolyBiofilm.cpp
  )

target_include_directories(Biofilm PUBLIC includes/Biofilm)
target_link_libraries(Biofilm PRIVATE MathUtility Bacteria Grid IO Forces)
