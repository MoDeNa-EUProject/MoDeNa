cmake_minimum_required (VERSION 2.8)
project (initMoments C CXX Fortran)

if( CMAKE_VERSION VERSION_GREATER "3.0" )
  cmake_policy(SET CMP0042 OLD)
  cmake_policy(SET CMP0026 OLD)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(CMAKE_CXX_FLAGS "-Wall -g")
set(CMAKE_BUILD_TYPE Debug)
message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
find_package(LAPACK REQUIRED)

add_executable(initMoments initMoments.C)
target_link_libraries(initMoments ${LAPACK_LIBRARIES})
