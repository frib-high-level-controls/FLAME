include_directories(
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

add_executable(example_linear
  example_linear.cpp
)
target_link_libraries(example_linear
  uscsi_core
)
