include_directories(
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../src
  ${CMAKE_CURRENT_BINARY_DIR}
  ${HDF5_CXX_INCLUDE_DIR}
)

add_executable(example_sim
  sim.cpp
)
target_link_libraries(example_sim
  uscsi_core
)

add_executable(example_customsim
  customsim.cpp
)
target_link_libraries(example_customsim
  uscsi_core
)

# example binaries should not be installed
