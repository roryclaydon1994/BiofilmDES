cmake_minimum_required (VERSION 3.16)

project (BiofilmDES LANGUAGES CXX VERSION 1.0)

set(CMAKE_BUILD_TYPE Debug) # Change to Release to disable assertions

### Require out-of-source builds i.e. you need to build this in a separate build dir
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

# All possible simulation MACROS listed below - see README.md
# -D CHAINING -D DEBUG_CHAINING -D MOVE_3D -D RANDOM_SEED

# Change CMAKE_BUILD_TYPE to Debug and turn on the following for development
# set(DEBUGGING_FLAGS "-Wall -pedantic -g -pg")

set(SIM_DEFS "-D CHAINING -D RANDOM_SEED")
string(APPEND SIM_DEFS  " -fopenmp ")

# Combine the SIM_DEFs and your debugging flags
set(ADDITIONAL_FLAGS ${SIM_DEFS} ${DEBUGGING_FLAGS})
message(STATUS "FLAGS: ${ADDITIONAL_FLAGS}")

set(CMAKE_CXX_FLAGS
    ${CMAKE_CXX_FLAGS}
    ${ADDITIONAL_FLAGS}
    )

# Cmake creates ; separated lists, so need to replace the ; with ' '
string(REPLACE ";" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
message("Running with CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}" )

include_directories (
  includes
  Grid/includes/Grid
  )

message("From top level:")
message("  Source (general): ${PROJECT_SOURCE_DIR}")
message("  Source (BiofilmDES): ${BiofilmDES_SOURCE_DIR}")

add_subdirectory(includes)
add_subdirectory(IO)
add_subdirectory(MathUtility)
add_subdirectory(Bacteria)
add_subdirectory(Grid)
add_subdirectory(Forces)
add_subdirectory(Biofilm)
add_subdirectory(Main)

add_subdirectory(Tests)
