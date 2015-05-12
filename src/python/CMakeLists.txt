
find_program(PYTHON "python" REQUIRED)

option(WITH_PYTHON_SWIG "Build Library with SWIG support for Python option" OFF)
if(WITH_PYTHON_SWIG)
    find_package(SWIG REQUIRED)
    include(${SWIG_USE_FILE})

    set(SWIG_MODULE_pymodena_EXTRA_DEPS ${MODENA_SRC_DIR}/modena.h)

    swig_add_module(pymodena python pymodena.i)
    swig_link_libraries(pymodena modena ${PYTHON_LIBRARIES})

    get_property(PYMODENA_LIB_LOCATION TARGET _pymodena PROPERTY LOCATION)
endif()


set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
set(SETUP_PY    ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
file(GLOB DEPS "*.py")
set(OUTPUT      ${CMAKE_CURRENT_BINARY_DIR}/build/timestamp)

configure_file(${SETUP_PY_IN} ${SETUP_PY})

add_custom_command(OUTPUT ${OUTPUT}
                   COMMAND ${PYTHON} ${SETUP_PY} build
                   COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                   DEPENDS ${DEPS}
                  )

if(WITH_PYTHON_SWIG)
    add_custom_target(python ALL DEPENDS ${OUTPUT} _pymodena)
else()
    add_custom_target(python ALL DEPENDS ${OUTPUT})
endif()

install(
        CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install --prefix=${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY python)"
       )

