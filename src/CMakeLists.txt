include_directories(
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${HDF5_CXX_INCLUDE_DIR}
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
  scsi/util.h
  scsi/base.h
  scsi/config.h
  scsi/constants.h

  scsi/state/vector.h
  scsi/state/matrix.h
  scsi/linear.h
  scsi/moment.h
  scsi/moment2.h
  scsi/moment2_sup.h
  scsi/cav_ttf.h
  scsi/rf_cavity.h
  scsi/chg_stripper.h

  scsi/h5loader.h
  scsi/h5writer.h
)
install(FILES ${inst_HEADERS}
  DESTINATION include/uscsi
)

add_library(uscsi_core SHARED
  base.cpp
  config.cpp

  linear.cpp

  moment.cpp

  moment2.cpp

  moment2_sup.cpp

  cav_ttf.cpp

  rf_cavity.cpp

  chg_stripper.cpp

  glps_parser.cpp glps_parser.h
  glps_ops.cpp
  glps.par.c glps.par.h
  glps.tab.c glps.tab.h

  h5loader.cpp
  h5writer.cpp

  ${inst_HEADERS}
)
target_link_libraries(uscsi_core
  ${Boost_THREAD_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${HDF5_LIBRARIES}
)
set_target_properties(uscsi_core
  PROPERTIES
    SOVERSION 0
)

install(TARGETS uscsi_core
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)


add_executable(test_lex
  test_lex.cpp
)
target_link_libraries(test_lex
  uscsi_core
)

add_executable(test_jb_1
  test_jb_1.cpp
)
target_link_libraries(test_jb_1
  uscsi_core
)

add_executable(test_jb_2
  test_jb_2.cpp
)
target_link_libraries(test_jb_2
  uscsi_core
)

add_executable(test_parse
  test_parse.cpp
)
target_link_libraries(test_parse
  uscsi_core
)

add_executable(test_config
  test_config.cpp
)
add_test(config test_config)
target_link_libraries(test_config
  uscsi_core
  ${Boost_PRG_EXEC_MONITOR_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

add_executable(h5_loader
  h5loadertest.cpp
)
target_link_libraries(h5_loader
  uscsi_core
)

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
