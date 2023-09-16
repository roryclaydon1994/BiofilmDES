cmake_minimum_required (VERSION 3.16)

add_library(IO STATIC includes/IO source/IO.cpp)
target_include_directories(IO PUBLIC includes/IO ${CMAKE_BINARY_DIR}/includes)
target_link_libraries(IO PRIVATE MathUtility Bacteria)
