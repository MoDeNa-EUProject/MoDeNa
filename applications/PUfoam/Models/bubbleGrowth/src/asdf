cmake_minimum_required (VERSION 2.6)
project (tutorialModels C CXX Fortran)

if( CMAKE_VERSION VERSION_GREATER "3.0" )
  cmake_policy(SET CMP0042 OLD)
  cmake_policy(SET CMP0026 OLD)
endif()

find_package(MODENA REQUIRED)

include_directories(${MODENA_INCLUDE_DIRS})
link_directories(${MODENA_LIBRARY_DIRS})

add_subdirectory(bubbleGrowth)
add_subdirectory(CFD_tool/0D)
