include_directories(
  ${PYTHON_INCLUDE_PATH}
  ${NUMPY_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../../src
)

add_library(_internal SHARED
  modmain.cpp
  modconfig.cpp
  modmachine.cpp
  modstate.cpp
)

target_link_libraries(_internal
  uscsi_core
  ${Boost_PYTHON_LIBRARY}
)

set_target_properties(_internal
  PROPERTIES
    PREFIX ""
)
