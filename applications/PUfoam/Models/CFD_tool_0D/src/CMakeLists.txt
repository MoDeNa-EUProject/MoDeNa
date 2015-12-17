cmake_minimum_required (VERSION 2.8)
project (QmomKinetics C CXX Fortran)

if( CMAKE_VERSION VERSION_GREATER "3.0" )
  cmake_policy(SET CMP0042 OLD)
  cmake_policy(SET CMP0026 OLD)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(CMAKE_CXX_FLAGS "-Wall -g")
set(CMAKE_BUILD_TYPE Debug)
message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
add_subdirectory(eigen)

#set(Boost_DEBUG 1)
find_package(Boost REQUIRED)
add_definitions( "-DHAS_BOOST" )
include_directories(${Boost_INCLUDE_DIRS})

find_package(MODENA REQUIRED)
include_directories(${MODENA_INCLUDE_DIRS})

add_executable(QmomKinetics QmomKinetics.cpp)
target_link_libraries(QmomKinetics ${Boost_LIBRARIES} eigen MODENA::modena)
