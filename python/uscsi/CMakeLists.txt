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
  ioc.py
)
set(PY_SRC_TEST
  test/__init__.py
  test/test_linear.py
  test/test_moment.py
  test/test_config.py
  test/test_jb.py
  test/moment_jb.lat
  test/moment_jb_2.lat
  test/parse1.lat
  test/latticeout_IMP_withPV_consolidate.lat
  test/latticeout_IMP_withPV_consolidate2.lat
  test/LS1.lat
  test/to_strl.lat
  test/to_strl_err.lat
  test/to_chg_str_err.lat
  test/latticeout_IMP_All.lat
  test/test_h5.lat
  test/test.h5
  test/jb.py
  test/jb_example.py
  test/jb_chg_str.py
  test/FrontEnd.lat

  test/data/axisData_41.txt
  test/data/axisData_85.txt
  test/data/Multipole41/thinlenlon_41.txt
  test/data/Multipole85/thinlenlon_85.txt
  test/data/Multipole41/CaviMlp_41.txt
  test/data/Multipole85/CaviMlp_85.txt
)

set(PY_SRC_DATA
  data/st.cmd
  data/core.db
  data/orbit.db
  data/common.template
  data/source.template
  data/quadrupole.template
  data/solenoid.template
)

install(FILES ${PY_SRC}
  DESTINATION ${Python_MODULE_DIR}/uscsi
)
install(FILES ${PY_SRC_TEST}
  DESTINATION ${Python_MODULE_DIR}/uscsi/test
)
install(FILES ${PY_SRC_DATA}
  DESTINATION ${Python_MODULE_DIR}/uscsi/data
)

# Copy .py files to the binary dir in order to run unit-tests

foreach(pyfile IN LISTS PY_SRC PY_SRC_TEST PY_SRC_DATA)
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
