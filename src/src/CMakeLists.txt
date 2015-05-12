# Build C library
add_library(modena SHARED modena.c)
target_link_libraries(modena ${LTDL_LIBRARIES} ${PYTHON_LIBRARIES})

#include(GenerateExportHeader)
#generate_export_header(modena)
set_property(TARGET modena PROPERTY VERSION ${PACKAGE_VERSION})
set_property(TARGET modena PROPERTY SOVERSION 1)
set_property(TARGET modena PROPERTY INTERFACE_MODENA_MAJOR_VERSION 1)
set_property(
  TARGET modena APPEND PROPERTY
  COMPATIBLE_INTERFACE_STRING MODENA_MAJOR_VERSION
)
set_target_properties(modena PROPERTIES PUBLIC_HEADER modena.h)

# Build Fortran library
add_library(fmodena SHARED fmodena.f90)
target_link_libraries(fmodena ${LTDL_LIBRARIES} ${PYTHON_LIBRARIES})

# Build executables
add_executable(twoTanksFullProblem twoTanksFullProblem.C)

# Build configuration files
configure_file(libmodena.pc.in libmodena.pc @ONLY)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/MODENAConfigVersion.cmake
  VERSION ${PACKAGE_VERSION}
  COMPATIBILITY AnyNewerVersion
)

export(
  TARGETS twoTanksFullProblem modena fmodena
  FILE ${CMAKE_CURRENT_BINARY_DIR}/MODENATargets.cmake
  NAMESPACE MODENA::
)

configure_file(
  cmake/MODENAConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/MODENAConfig.cmake
  COPYONLY
)

# Installation
install(
  TARGETS twoTanksFullProblem modena fmodena
  EXPORT MODENATargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static
  INCLUDES DESTINATION include
  PUBLIC_HEADER DESTINATION include
)

install(FILES libmodena.pc DESTINATION lib/pkgconfig)
install(FILES fmodena.mod DESTINATION include)

install(
  EXPORT MODENATargets
  FILE MODENATargets.cmake
  NAMESPACE MODENA::
  DESTINATION ${ConfigPackageLocation}
)

install(
  FILES cmake/MODENAConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/MODENAConfigVersion.cmake
  DESTINATION ${ConfigPackageLocation}
)

