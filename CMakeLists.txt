project("ultrascsi" C CXX)

cmake_minimum_required(VERSION 2.8)

option(WERROR "Treat warnings as errors")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE DEBUG) # default to debug build
endif()
if(CMAKE_COMPILER_IS_GNUCC)
  # Turn on extra GCC warnings
  add_definitions(-Wall -Wextra -Wno-unused-parameter)
  if(WERROR)
    add_definitions(-Werror)
  endif()
endif()

find_package(Boost 1.55.0 REQUIRED)

add_subdirectory(src)
add_subdirectory(test)
