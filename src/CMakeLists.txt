cmake_minimum_required (VERSION 3.0)
project (modena VERSION 1.0 LANGUAGES C CXX Fortran)

message("*************************************************")
message("-- MoDeNa Version: ${PROJECT_VERSION}")
message("-- MoDeNa Build Dir: ${CMAKE_CURRENT_BINARY_DIR}")
message("-- MoDeNa Install Prefix: ${CMAKE_INSTALL_PREFIX}")
message("*************************************************")


# Outlaw "in-source" build -------------------------------------------------- #
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)


set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)

# Build Types --------------------------------------------------------------- #

set(CMAKE_BUILD_TYPE Debug)

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

# Standardise Install Locations --------------------------------------------- #
set(INSTALL_LIB_DIR     lib     CACHE PATH "Library install directories")
set(INSTALL_BIN_DIR     bin     CACHE PATH "Executables install directories")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Header files install directory")

foreach(p LIB BIN INCLUDE CMAKE)                 # Make relative paths absolute
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

# --------------------------------------------------------------------------- #

set(MODENA_SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/src)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

# Locate Libraries ---------------------------------------------------------- #

include(CheckCInline)
check_c_inline(C_INLINE)

find_package(LTDL REQUIRED)

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})


include_directories(src)

set(ConfigPackageLocation lib/cmake/modena)
install(DIRECTORY cmake/ DESTINATION ${ConfigPackageLocation})

add_subdirectory(src)
add_subdirectory(python)

