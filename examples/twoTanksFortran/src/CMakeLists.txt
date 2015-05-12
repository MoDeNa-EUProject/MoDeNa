cmake_minimum_required (VERSION 2.6)
project (tutorialModels C CXX Fortran)

if( CMAKE_VERSION VERSION_GREATER "3.0" )
  cmake_policy(SET CMP0042 OLD)
  cmake_policy(SET CMP0026 OLD)
endif()



set(CMAKE_BUILD_TYPE Release)

find_package(MODENA REQUIRED)

add_executable(flowRateExact flowRateExact.C)

add_executable(twoTanksMacroscopicProblemFortran twoTanksMacroscopicProblemFortran.f90)
target_link_libraries(twoTanksMacroscopicProblemFortran MODENA::modena)

