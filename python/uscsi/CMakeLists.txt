include_directories(
  ${Boost_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_PATH}
  ${NUMPY_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../../src
)

add_library(_internal SHARED
  modmain.cpp
  modconfig.cpp
  modmachine.cpp
  modstate.cpp
  pyscsi.h
)

target_link_libraries(_internal
  uscsi_core
  ${Boost_PYTHON_LIBRARY}
)

set_target_properties(_internal
  PROPERTIES
    PREFIX ""
)

set(PY_SRC
  __init__.py
  test/__init__.py
  test/test_linear.py
  test/test_config.py
)

# Copy .py files to the binary dir in order to run unit-tests

foreach(pyfile IN LISTS PY_SRC)
  add_custom_command(
    OUTPUT ${pyfile}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${pyfile} ${pyfile}
    MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${pyfile}
    DEPENDS ${_pydir}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  list(APPEND PY_BIN ${pyfile})
endforeach()

add_custom_target(testdata ALL
  DEPENDS ${PY_BIN}
)

if(HAVE_NOSE)
  add_test(
    NAME nosetest
    COMMAND ${PYTHON_EXECUTABLE} -m nose.core uscsi
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
  )
endif()
