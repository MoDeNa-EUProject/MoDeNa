cmake_minimum_required (VERSION 2.8)
project (tutorialModels C Fortran)

set (CMAKE_NO_SYSTEM_FROM_IMPORTED yes)

if( CMAKE_VERSION VERSION_GREATER "3.0" )
  cmake_policy(SET CMP0042 OLD)
  cmake_policy(SET CMP0026 OLD)
endif()

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
find_library(FSON fson PATHS $ENV{HOME}/lib)

include_directories(${MODENA_INCLUDE_DIRS})
include_directories($ENV{HOME}/include)
link_directories(${MODENA_LIBRARY_DIRS})

set (CMAKE_Fortran_FLAGS "-ffree-line-length-none -O3")
set (CMAKE_Fortran_MODULE_DIRECTORY mod)
file (GLOB _sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} src/*.f*)
add_executable(degas ${_sources})
target_link_libraries(degas MODENA::modena ${FSON})
