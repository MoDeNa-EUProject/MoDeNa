find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

set (CMAKE_Fortran_FLAGS "-ffree-line-length-none -O3")
set (CMAKE_Fortran_MODULE_DIRECTORY mod)
file (GLOB _sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} src/*.f*)
add_executable(kfoam ${_sources})
target_link_libraries(kfoam ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MODENA::modena)
