include_directories(
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

add_executable(run_uscsi
  cli.cpp
)
target_link_libraries(run_uscsi
  uscsi_core
)
