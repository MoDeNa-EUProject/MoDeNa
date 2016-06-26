cmake_minimum_required (VERSION 2.8)
project (bubbleGrowth C CXX Fortran)

if( CMAKE_VERSION VERSION_GREATER "3.0" )
    cmake_policy(SET CMP0042 OLD)
    cmake_policy(SET CMP0026 OLD)
endif()

find_package(MODENA QUIET)

set (CMAKE_NO_SYSTEM_FROM_IMPORTED yes)

set (CMAKE_Fortran_FLAGS "-ffree-line-length-none -O3")

set (CMAKE_Fortran_MODULE_DIRECTORY mod)
file (GLOB _sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} src/*.f*)
add_executable(bblgrExact ${_sources})

if(MODENA_FOUND)
    message("MoDeNa found")
    include_directories(${MODENA_INCLUDE_DIRS})
    link_directories(${MODENA_LIBRARY_DIRS})
    target_link_libraries(bblgrExact MODENA::modena)
else()
    message("MoDeNa not found")
    file (GLOB _fms RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} src/fake_modena.c)
    message(${_sources})
    message(${CMAKE_CURRENT_SOURCE_DIR})
    message(${_fms})
    add_library(fake_modena STATIC ${_fms})
    target_link_libraries(bblgrExact fake_modena)
endif()
