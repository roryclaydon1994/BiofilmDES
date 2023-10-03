set(OUTPUT_DIR "${PROJECT_SOURCE_DIR}/../GeneratedOutput/SimOutput/" )

message("OUTPUT_DIR: ${OUTPUT_DIR}")
set(SIM_DIR ${OUTPUT_DIR})
message("Old outpath: ${CMAKE_CURRENT_SOURCE_DIR}/outputPaths.hpp.in")
message("Move outpath: ${CMAKE_CURRENT_BINARY_DIR}/outputPaths.hpp")
configure_file(
  outputPaths.hpp.in
  outputPaths.hpp @ONLY
  )
message("include: ${CMAKE_CURRENT_BINARY_DIR}")
include_directories(${CMAKE_CURRENT_BINARY_DIR})

if(EXISTS ${OUTPUT_DIR})
  message("${OUTPUT_DIR} already exists")
else()
  message("Creating simulation outut directory at: ${OUTPUT_DIR}")
  file(MAKE_DIRECTORY ${OUTPUT_DIR})
endif()
