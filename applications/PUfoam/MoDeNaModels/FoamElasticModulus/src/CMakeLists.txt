set (CMAKE_C_FLAGS "-O3 -w")
file (GLOB _sources *.c*)
add_executable(spherepack ${_sources})
target_link_libraries(spherepack m)
