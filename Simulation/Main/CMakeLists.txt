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
  project (Main )
else()
  cmake_policy(SET CMP0048 NEW)
  project (Main LANGUAGES CXX VERSION 1.0)
endif()

# specify the C++ standard
#set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories (includes)
#configure_file(includes/Version.h.in includes/Version.h)

add_executable(main.out main.cpp)
target_link_libraries(
  main.out
  PRIVATE
  MathUtility
  Forces
  Phage
  Bacteria
  Biofilm
  IO
)

#add_executable(findOrientationalStiffness.out findOrientationalStiffness.cpp)
#target_compile_definitions( findOrientationalStiffness.out PUBLIC CHAINING RANDOM_SEED )
#target_link_libraries(
#  findOrientationalStiffness.out
#  PRIVATE
#  MathUtility
#  Forces
#  Phage
#  Bacteria
#  Biofilm
#  IO
#)

#add_executable(findCriticalChainingProb.out findCriticalChainingProb.cpp)
#target_compile_definitions( findCriticalChainingProb.out PUBLIC CHAINING RANDOM_SEED )
#target_link_libraries(
#  findCriticalChainingProb.out
#  PRIVATE
#  MathUtility
#  Forces
#  Phage
#  Bacteria
#  Biofilm
#  IO
#)

#add_executable(findLcDist.out findLcDist.cpp)
#target_compile_definitions( findLcDist.out PUBLIC CHAINING RANDOM_SEED )
#target_link_libraries(
#  findLcDist.out
#  PRIVATE
#  MathUtility
#  Forces
#  Phage
#  Bacteria
#  Biofilm
#  IO
#)

#add_executable(findScaling.out findScaling.cpp)
#target_compile_definitions( findScaling.out PUBLIC CHAINING RANDOM_SEED )
#target_link_libraries(
#  findScaling.out
#  PRIVATE
#  MathUtility
#  Forces
#  Phage
#  Bacteria
#  Biofilm
#  IO
#)

#add_executable(findOrientationalStiffness.out findOrientationalStiffness.cpp)
#target_compile_definitions( findOrientationalStiffness.out PUBLIC CHAINING RANDOM_SEED )
#target_link_libraries(
#  findOrientationalStiffness.out
#  PRIVATE
#  MathUtility
#  Forces
#  Phage
#  Bacteria
#  Biofilm
#  IO
#)
