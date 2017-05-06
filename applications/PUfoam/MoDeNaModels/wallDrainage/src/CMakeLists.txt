cmake_minimum_required(VERSION 2.8)
set (CMAKE_FORTRAN_FLAGS "-O3")
file (GLOB _sources *.f*)
set (CMAKE_Fortran_MODULE_DIRECTORY mod)

find_package(FSON REQUIRED)
include_directories(${FSON_INCLUDE_DIRS})

find_package(BSPLINE REQUIRED)
include_directories(${BSPLINE_INCLUDE_DIRS})

add_executable(wallDrainage ${_sources})
target_link_libraries(
    wallDrainage
    ${FSON_LIBRARIES}
    ${BSPLINE_LIBRARIES}
)
