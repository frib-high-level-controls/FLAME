include_directories(
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  # ${HDF5_CXX_INCLUDE_DIR}
)

add_custom_command(
  OUTPUT glps.par.c glps.par.h
  COMMAND ${YACC_EXECUTABLE}
  ARGS
    -o glps.par.c
    --defines=glps.par.h
    ${CMAKE_CURRENT_SOURCE_DIR}/glps.y
  MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/glps.y
)

add_custom_command(
  OUTPUT glps.tab.c glps.tab.h
  COMMAND ${LEX_EXECUTABLE}
  ARGS
    -o glps.tab.c
    --header-file=glps.tab.h
    ${CMAKE_CURRENT_SOURCE_DIR}/glps.l
  MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/glps.l
)

set_source_files_properties(
    glps.par.c glps.par.h
    glps.tab.c glps.tab.h
  PROPERTIES
    GENERATED TRUE
)

set(inst_HEADERS
  flame/util.h
  flame/base.h
  flame/config.h
  flame/constants.h

  flame/state/vector.h
  flame/state/matrix.h
  flame/linear.h
  flame/moment.h
  flame/moment_sup.h
  flame/rf_cavity.h
  flame/chg_stripper.h

  # flame/h5loader.h
  # flame/h5writer.h
)
install(FILES ${inst_HEADERS}
  DESTINATION include/flame
)

add_library(flame_core SHARED
  base.cpp
  config.cpp

  linear.cpp
  moment.cpp
  moment_sup.cpp
  rf_cavity.cpp
  chg_stripper.cpp

  glps_parser.cpp glps_parser.h
  glps_ops.cpp
  glps.par.c glps.par.h
  glps.tab.c glps.tab.h

  util.cpp

  # h5loader.cpp
  # h5writer.cpp

  ${inst_HEADERS}
)
target_link_libraries(flame_core
  ${Boost_THREAD_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  # ${HDF5_LIBRARIES}
)
set_target_properties(flame_core
  PROPERTIES
    SOVERSION 0
)

install(TARGETS flame_core
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)


add_executable(test_lex
  test_lex.cpp
)
target_link_libraries(test_lex
  flame_core
)

add_executable(test_jb_2
  test_jb_2.cpp
)
target_link_libraries(test_jb_2
  flame_core
)

add_executable(test_parse
  test_parse.cpp
)
target_link_libraries(test_parse
  flame_core
)

add_executable(test_config
  test_config.cpp
)
add_test(config test_config)
target_link_libraries(test_config
  flame_core
  ${Boost_PRG_EXEC_MONITOR_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

add_executable(test_util
  test_util.cpp
)
add_test(util test_util)
target_link_libraries(test_util
  flame_core
  ${Boost_PRG_EXEC_MONITOR_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

# add_executable(h5_loader
#   h5loadertest.cpp
# )
# target_link_libraries(h5_loader
#  flame_core
# )

if(UNIX)

  add_test(recurse1
    /bin/sh ${CMAKE_CURRENT_SOURCE_DIR}/recurse-test.sh
      ./test_parse
      ${CMAKE_CURRENT_SOURCE_DIR}/data/tracy_1.lat
      ${CMAKE_CURRENT_SOURCE_DIR}/data/tracy_1_out.lat
  )

  add_test(recurse2
    /bin/sh ${CMAKE_CURRENT_SOURCE_DIR}/recurse-test.sh
      ./test_parse
      ${CMAKE_CURRENT_SOURCE_DIR}/data/tracy_2.lat
      ${CMAKE_CURRENT_SOURCE_DIR}/data/tracy_2_out.lat
  )

endif()
