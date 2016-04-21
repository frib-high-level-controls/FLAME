include_directories(
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

add_executable(run_flame
  cli.cpp
)
target_link_libraries(run_flame
  uscsi_core
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
)

install(TARGETS run_flame
  RUNTIME DESTINATION bin
)
