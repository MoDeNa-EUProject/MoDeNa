cmake_minimum_required (VERSION 2.6)
project (tutorialModels C Fortran)

if( CMAKE_VERSION VERSION_GREATER "3.0" )
  cmake_policy(SET CMP0042 OLD)
  cmake_policy(SET CMP0026 OLD)
endif()

set(CMAKE_BUILD_TYPE Release)

add_custom_target(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

add_custom_target(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

find_package(MODENA REQUIRED)

#simple kinetics
set (CMAKE_Fortran_FLAGS "-ffree-line-length-none -O3")
set (CMAKE_Fortran_MODULE_DIRECTORY simpleKinetics_mod)
file (GLOB _simpleKinetics_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} simpleKinetics_src/*.f*)
add_executable(simpleKinetics ${_simpleKinetics_sources})
target_link_libraries(simpleKinetics MODENA::modena)
