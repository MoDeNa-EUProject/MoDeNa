set (CMAKE_FORTRAN_FLAGS "-O3")
file (GLOB _sources *.f*)
set (CMAKE_Fortran_MODULE_DIRECTORY mod)
add_executable(wallDrainage ${_sources})
