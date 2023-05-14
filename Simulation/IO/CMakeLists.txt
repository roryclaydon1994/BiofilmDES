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
  project (IO )
else()
  cmake_policy(SET CMP0048 NEW)
  project (IO LANGUAGES CXX VERSION 1.0)
endif()


# specify the C++ standard
#set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(IO STATIC includes/IO source/IO.cpp)
target_include_directories(IO PUBLIC includes/IO ${CMAKE_BINARY_DIR}/includes)
target_link_libraries(IO PRIVATE MathUtility Bacteria)
