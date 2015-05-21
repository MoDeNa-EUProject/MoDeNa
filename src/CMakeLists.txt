cmake_minimum_required (VERSION 2.6)
project (modena C CXX Fortran)

if( CMAKE_VERSION VERSION_GREATER "3.0" )
    cmake_policy(SET CMP0042 OLD)
    cmake_policy(SET CMP0026 OLD)
endif()

if(APPLE)
    set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
endif(APPLE)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)

set(PACKAGE_VERSION 0.5)
set(CMAKE_BUILD_TYPE Release)

set(MODENA_SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/src)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

find_package(LTDL REQUIRED)

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

include_directories(src)

set(ConfigPackageLocation lib/cmake/modena)
install(DIRECTORY cmake/ DESTINATION ${ConfigPackageLocation})

add_subdirectory(src)
add_subdirectory(python)

