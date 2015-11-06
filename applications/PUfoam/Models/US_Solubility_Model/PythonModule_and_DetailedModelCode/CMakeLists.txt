cmake_minimum_required (VERSION 2.6)
project (tutorialModels C Fortran)

if( CMAKE_VERSION VERSION_GREATER "3.0" )
  cmake_policy(SET CMP0042 OLD)
  cmake_policy(SET CMP0026 OLD)
endif()

set (CMAKE_Fortran_FLAGS "-fdefault-real-8")
set (CMAKE_Fortran_MODULE_DIRECTORY mod)
file (GLOB _sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DetailedModelCode/*.f*)
add_executable(PCSAFT_Henry ${_sources})
