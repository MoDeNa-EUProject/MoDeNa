cmake_minimum_required (VERSION 2.8)
project (modena C CXX Fortran)

if( CMAKE_VERSION VERSION_GREATER "3.0" )
    cmake_policy(SET CMP0042 OLD)
    cmake_policy(SET CMP0026 OLD)
    cmake_policy(SET CMP0028 OLD)
endif()


set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)

set(PACKAGE_VERSION 1.0)

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

set(MODENA_SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/src)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

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

