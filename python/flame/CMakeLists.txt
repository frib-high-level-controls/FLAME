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
  pyflame.h
)
install(TARGETS _internal
  LIBRARY DESTINATION ${Python_MODULE_DIR}/flame
)

target_link_libraries(_internal
  flame_core
)

set_target_properties(_internal
  PROPERTIES
    PREFIX ""
    NO_SONAME YES
)

set(PY_SRC
  __init__.py
  __main__.py
  ioc.py

  test/__init__.py
  test/test_linear.py
  test/test_moment.py
  test/test_config.py
  test/test_beam_dynamics.py
)
set(PY_DATA
  # data for test_config.py
  test/test_h5.lat
  test/test.h5

  # lattice files for unit tests
  test/parse1.lat
  test/LS1.lat
  test/to_strl.lat
  test/to_chg_str_err.lat
  test/FE_latticeE.lat
  test/LS1FS1_latticeE.lat
  test/ALL_latticeE.lat
  test/Arc_Ds.lat
  test/LS1_Gen.lat
  test/LS2_Gen.lat
  test/Test7GapCav.lat

  # cavity data for test lattice files
  test/data/axisData_29.txt
  test/data/axisData_41.txt
  test/data/axisData_53.txt
  test/data/axisData_85.txt
  test/data/thinlenlon_29.lat
  test/data/thinlenlon_41.lat
  test/data/thinlenlon_53.lat
  test/data/thinlenlon_85.lat
  test/data/thinlenlon_7g.lat
  test/data/Multipole29/thinlenlon_29.txt
  test/data/Multipole41/thinlenlon_41.txt
  test/data/Multipole53/thinlenlon_53.txt
  test/data/Multipole85/thinlenlon_85.txt
  test/data/Multipole29/CaviMlp_29.txt
  test/data/Multipole41/CaviMlp_41.txt
  test/data/Multipole53/CaviMlp_53.txt
  test/data/Multipole85/CaviMlp_85.txt
  

  # vaccel IOC template
  data/st.cmd
  data/core.db
  data/orbit.db
  data/common.template
  data/source.template
  data/drift.template
  data/rfcavity.template
  data/quadrupole.template
  data/sextupole.template
  data/sbend.template
  data/solenoid.template
  data/edipole.template
  data/equad.template
  data/stripper.template
  data/bpm.template
  data/marker.template
  data/orbtrim.template
)

# Copy .py files to the binary dir in order to run unit-tests

foreach(pyfile IN LISTS PY_SRC PY_DATA)
  get_filename_component(_dir ${pyfile} DIRECTORY)
  add_custom_command(
    OUTPUT ${pyfile}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${pyfile} ${pyfile}
    MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${pyfile}
    DEPENDS ${_pydir}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  install(FILES ${pyfile}
    DESTINATION ${Python_MODULE_DIR}/flame/${_dir}
  )
endforeach()

add_custom_target(pyfiles ALL
  DEPENDS ${PY_SRC} ${PY_DATA}
)

if(HAVE_NOSE)
  add_test(
    NAME nosetest
    COMMAND ${PYTHON_EXECUTABLE} -m nose.core flame
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
  )
endif()
