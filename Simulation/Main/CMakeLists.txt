cmake_minimum_required (VERSION 3.16)

include_directories (includes)

add_executable(main.out main.cpp)

target_link_libraries(
  main.out
  PRIVATE
  MathUtility
  Forces
  Bacteria
  Biofilm
  IO
)