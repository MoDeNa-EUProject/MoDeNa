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

# Set search path for SUNDIALS libraries
SET(SUNDIALS_LIB_DIR /usr/local/lib)

# Find the SUNDIALS solver's library
FIND_LIBRARY(SUNDIALS_SOLVER_LIB
  sundials_cvode ${SUNDIALS_LIB_DIR}
  DOC "CVODE library")

FIND_LIBRARY(SUNDIALS_SOLVER_FLIB
  sundials_fcvode ${SUNDIALS_LIB_DIR}
  DOC "CVODE F77-C library")

# Find the NVECTOR library
FIND_LIBRARY(SUNDIALS_NVEC_LIB
  sundials_nvecserial ${SUNDIALS_LIB_DIR}
  DOC "NVECTOR library")

FIND_LIBRARY(SUNDIALS_NVEC_FLIB
  sundials_fnvecserial ${SUNDIALS_LIB_DIR}
  DOC "NVECTOR F77-C library")

# Set additional libraries
SET(SUNDIALS_EXTRA_LIB  -lm /usr/lib/x86_64-linux-gnu/librt.so
    CACHE STRING "Additional libraries")

# Set Blas/Lapack libraries
# SET(LAPACK_LIBRARIES  CACHE STRING "Lapack libraries")
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

# List of all libraries
SET(SUNDIALS_LIBS
  ${SUNDIALS_SOLVER_FLIB} ${SUNDIALS_SOLVER_LIB}
  ${SUNDIALS_NVEC_FLIB} ${SUNDIALS_NVEC_LIB}
  ${SUNDIALS_EXTRA_LIB})
IF(LAPACK_FOUND)
  LIST(APPEND SUNDIALS_LIBS ${LAPACK_LIBRARIES})
ENDIF(LAPACK_FOUND)

include_directories($ENV{HOME}/include)
find_library(FSON fson PATHS $ENV{HOME}/lib)
find_library(BSPLINE bspline PATHS $ENV{HOME}/lib)
if(MODENA_FOUND)
    message("MoDeNa found")
    include_directories(${MODENA_INCLUDE_DIRS})
    link_directories(${MODENA_LIBRARY_DIRS})
    target_link_libraries(
        bblgrExact MODENA::modena ${FSON} ${BSPLINE} ${SUNDIALS_LIBS})
else()
    message("MoDeNa not found")
    file (GLOB _fms RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} src/fake_modena.c)
    message(${_sources})
    message(${CMAKE_CURRENT_SOURCE_DIR})
    message(${_fms})
    add_library(fake_modena STATIC ${_fms})
    target_link_libraries(
        bblgrExact fake_modena ${FSON} ${BSPLINE} ${SUNDIALS_LIBS})
endif()
