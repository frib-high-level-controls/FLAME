include_directories(
  ${Boost_INCLUDE_DIRS}
  ${Python_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../../src
)

add_library(_internal SHARED
  modmain.cpp
  modconfig.cpp
  modmachine.cpp
  modstate.cpp
  pyscsi.h
)
install(TARGETS _internal
  LIBRARY DESTINATION ${Python_MODULE_DIR}/uscsi
)

target_link_libraries(_internal
  uscsi_core
)

set_target_properties(_internal
  PROPERTIES
    PREFIX ""
    NO_SONAME YES
)

set(PY_SRC
  __init__.py
)
set(PY_SRC_TEST
  test/__init__.py
  test/test_linear.py
  test/test_moment.py
  test/test_config.py
  test/test_jb.py
  test/moment_jb.lat
  test/moment_jb_2.lat
  test/latticeout_IMP_withPV_consolidate.lat
  test/test_h5.lat
  test/test.h5
)

install(FILES ${PY_SRC}
  DESTINATION ${Python_MODULE_DIR}/uscsi
)
install(FILES ${PY_SRC_TEST}
  DESTINATION ${Python_MODULE_DIR}/uscsi/test
)

# Copy .py files to the binary dir in order to run unit-tests

foreach(pyfile IN LISTS PY_SRC PY_SRC_TEST)
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
