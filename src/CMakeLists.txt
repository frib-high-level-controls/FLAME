include_directories(
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
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

add_library(uscsi_core SHARED
  scsi/util.h
  scsi/base.h
  base.cpp
  config.cpp scsi/config.h

  scsi/state/vector.h
  scsi/state/matrix.h

  linear.cpp
  scsi/linear.h

  moment.cpp
  scsi/moment.h

  glps_parser.cpp glps_parser.h
  glps_ops.cpp
  glps.par.c glps.par.h
  glps.tab.c glps.tab.h
)

add_executable(test_lex
  test_lex.cpp
)
target_link_libraries(test_lex
  uscsi_core
)

add_executable(test_parse
  test_parse.cpp
)
target_link_libraries(test_parse
  uscsi_core
)
