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
project (Particle LANGUAGES CXX VERSION 1.0)

# specify the C++ standard
#set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(Particle STATIC source/particle.cpp includes/Particle/particle.hpp)
target_include_directories(Particle PUBLIC includes/Particle)
target_link_libraries(Particle PRIVATE MathUtility)
