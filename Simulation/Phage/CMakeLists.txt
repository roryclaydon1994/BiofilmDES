cmake_minimum_required (VERSION 3.16)

add_library(Phage STATIC source/Phage.cpp includes/Phage/Phage.hpp)

target_include_directories(Phage PUBLIC includes/Phage)

target_link_libraries(Phage PRIVATE MathUtility Bacteria Grid)
