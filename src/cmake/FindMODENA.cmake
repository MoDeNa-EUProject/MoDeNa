# Read-Only variables:
#  MODENA_FOUND - system has the MODENA library
#  MODENA_INCLUDE_DIR - the MODENA include directory
#  MODENA_LIBRARIES - The libraries needed to use MODENA
#  MODENA_VERSION - This is set to $major.$minor.$revision$path (eg. 0.4.1)

if (UNIX)
  find_package(PkgConfig QUIET)
  pkg_check_modules(_MODENA QUIET libmodena)
  set(MODENA_DEFINITIONS ${_MODENA_CFLAGS_OTHER})
endif ()

find_path(MODENA_INCLUDE_DIR
  NAMES
    modena.h
  HINTS
    ${MODENA_ROOT_DIR}
    ${_MODENA_INCLUDEDIR}
  PATH_SUFFIXES
    include
)

find_package(LTDL)
find_package(PythonLibs)

set(MODENA_INCLUDE_DIR "${MODENA_INCLUDE_DIR}")
set(MODENA_INCLUDE_DIRS ${MODENA_INCLUDE_DIR} ${LTDL_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH})

if(WIN32 AND NOT CYGWIN)
  if(MSVC)
    find_library(MODENA
      NAMES
        "modena"
      HINTS
        ${MODENA_ROOT_DIR}
      PATH_SUFFIXES
        bin
        lib
    )

    mark_as_advanced(MODENA)
    set(MODENA_LIBRARIES ${MODENA} ws2_32)
  else()
      # bother supporting this?
  endif()
else()

  find_library(MODENA_LIBRARY
    NAMES
      modena libmodena
    HINTS
      ${_MODENA_LIBDIR}
    PATH_SUFFIXES
      lib
  )

  mark_as_advanced(MODENA_LIBRARY)

  set(MODENA_LIBRARIES ${MODENA_LIBRARY} ${LTDL_LIBRARIES} ${PYTHON_LIBRARIES})

endif()

if (MODENA_INCLUDE_DIR)
  if (_MODENA_VERSION)
     set(MODENA_VERSION "${_MODENA_VERSION}")
  elseif(MODENA_INCLUDE_DIR AND EXISTS "${MODENA_INCLUDE_DIR}/modena-version.h")
     file(STRINGS "${MODENA_INCLUDE_DIR}/modena-version.h" modena_version_str
        REGEX "^#define[\t ]+MODENA_VERSION[\t ]+\([0-9.]+\)[\t ]+$")

     string(REGEX REPLACE "^.*MODENA_VERSION[\t ]+\([0-9.]+\)[\t ]+$"
        "\\1" MODENA_VERSION "${modena_version_str}")
  endif ()
endif ()

include(FindPackageHandleStandardArgs)

if (MODENA_VERSION)
   find_package_handle_standard_args(MODENA
    REQUIRED_VARS
      MODENA_LIBRARIES
      MODENA_INCLUDE_DIR
    VERSION_VAR
      MODENA_VERSION
    FAIL_MESSAGE
      "Could NOT find MODENA version"
  )
else ()
   find_package_handle_standard_args(MODENA "Could NOT find MODENA uuuurh"
      MODENA_LIBRARIES
      MODENA_INCLUDE_DIR
  )
endif ()

list(REMOVE_DUPLICATES MODENA_LIBRARIES)
list(REMOVE_DUPLICATES MODENA_INCLUDE_DIRS)

mark_as_advanced(MODENA_INCLUDE_DIR MODENA_LIBRARIES)
