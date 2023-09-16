cmake_minimum_required (VERSION 3.16)

add_library(
  Forces STATIC
  source/forces.cpp
  includes/Forces/forces.hpp
  )
target_include_directories(Forces PUBLIC includes/Forces)
target_link_libraries(Forces PRIVATE MathUtility Bacteria Grid IO)