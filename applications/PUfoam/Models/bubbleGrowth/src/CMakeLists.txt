cmake_minimum_required (VERSION 2.8)
project (bubbleGrowth C CXX Fortran)

if( CMAKE_VERSION VERSION_GREATER "3.0" )
    cmake_policy(SET CMP0042 OLD)
    cmake_policy(SET CMP0026 OLD)
endif()

find_package(MODENA REQUIRED)

include_directories(${MODENA_INCLUDE_DIRS})
link_directories(${MODENA_LIBRARY_DIRS})

set (CMAKE_NO_SYSTEM_FROM_IMPORTED yes)

set (CMAKE_Fortran_FLAGS "-ffree-line-length-none -O3")

set (CMAKE_Fortran_MODULE_DIRECTORY mod)
file (GLOB _sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} src/*.f*)
add_executable(bblgrExact ${_sources})
target_link_libraries(bblgrExact MODENA::modena)
