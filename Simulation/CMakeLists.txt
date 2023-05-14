cmake_minimum_required (VERSION 3.16)
set(CMAKE_BUILD_TYPE Release)
message(" CMAKE_BUILD_TYPE IS ${CMAKE_BUILD_TYPE}")

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR
    "You cannot build in a source directory (or any directory with a CMakeLists.txt file).
    Please make a build subdirectory.
    Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# use -DNDEBUG to turn off assertions (speeds up runtime)
#SIM_DEFS="-D CHAINING -D DEBUG_CHAINING -D PHAGE -D MOVE_3D -D RANDOM_SEED"
#set(SIM_DEFS "-D AG43 -D RANDOM_SEED")
#set(DEBUGGING_FLAGS "-Wall -pedantic -g -pg")

set(SIM_DEFS "-D CHAINING -D RANDOM_SEED")
string(APPEND SIM_DEFS  " -fopenmp ")

set(ADDITIONAL_FLAGS ${SIM_DEFS} ${DEBUGGING_FLAGS})
message(STATUS "FLAGS: ${ADDITIONAL_FLAGS}")

project (BiofilmDES LANGUAGES CXX VERSION 1.0)
set(CMAKE_CXX_FLAGS
    ${CMAKE_CXX_FLAGS}
    ${ADDITIONAL_FLAGS}
    )
string(REPLACE ";" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
message("Running with CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}" )

#find_package(Boost 1.71.0 REQUIRED system)
include_directories (
  includes
  Grid/includes/Grid
  #"/usr/local/cuda-11.0/targets/x86_64-linux/include"
  )

message("From top level:")
message("  Source (general): ${PROJECT_SOURCE_DIR}")
message("  Source (BiofilmDES): ${BiofilmDES_SOURCE_DIR}")

add_subdirectory(includes)
add_subdirectory(IO)
add_subdirectory(MathUtility)
#add_subdirectory(Particle)
add_subdirectory(Phage)
add_subdirectory(Bacteria)
add_subdirectory(Grid)
add_subdirectory(Forces)
add_subdirectory(Biofilm)
add_subdirectory(Main)

#enable_testing()
#list(APPEND CMAKE_CTEST_ARGUMENTS "--verbose;--stop-on-failure;--output-on-failure")
#add_subdirectory(Tests)
#add_test(NAME TestBacteria COMMAND TestBac )
#add_test(NAME TestGrid COMMAND TestGrid )
#add_test(NAME TestPolyGrid COMMAND TestPolyGrid )
#add_test(NAME TestInteractions COMMAND TestInteractions )
#add_test(NAME TestPolyInteractions COMMAND TestPolyInteractions )
add_test(NAME TestPolyBiofilm COMMAND TestPolyBiofilm )
#add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} ${CMAKE_CTEST_ARGUMENTS})
