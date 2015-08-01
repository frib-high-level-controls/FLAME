include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

add_executable(test_linear
  test_linear.cpp
)
target_link_libraries(test_linear
  uscsi_core
)
