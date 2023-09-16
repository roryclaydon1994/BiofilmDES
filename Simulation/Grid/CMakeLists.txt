cmake_minimum_required (VERSION 3.16)

add_library(
  Grid
  INTERFACE
)
target_include_directories(Grid INTERFACE includes/Grid)
