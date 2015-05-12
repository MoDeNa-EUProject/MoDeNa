cmake_minimum_required (VERSION 2.6)
project (tutorialModels C CXX)

if( CMAKE_VERSION VERSION_GREATER "3.0" )
  cmake_policy(SET CMP0042 OLD)
  cmake_policy(SET CMP0026 OLD)
endif()


set(CMAKE_BUILD_TYPE Release)

find_package(MODENA REQUIRED)

add_executable(flowRateExact flowRateExact.C)
target_link_libraries(flowRateExact MODENA::modena)

add_executable(twoTanksMacroscopicProblem twoTanksMacroscopicProblem.C)
target_link_libraries(twoTanksMacroscopicProblem MODENA::modena)

