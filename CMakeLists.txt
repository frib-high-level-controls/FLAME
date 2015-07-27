cmake_minimum_required(VERSION 2.8.8)

project("uscsi" C CXX)

option(USE_WERROR "Treat warnings as errors")
option(USE_PYTHON "Build python support" ON)

message(STATUS "Option PYTHON ${USE_PYTHON}")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE DEBUG) # default to debug build
endif()
if(CMAKE_COMPILER_IS_GNUCC)
  # Turn on extra GCC warnings
  add_definitions(-Wall -Wextra -Wno-unused-parameter)
  if(USE_WERROR)
    add_definitions(-Werror)
  endif()
endif()

find_package(Boost 1.55.0 COMPONENTS system)

if(NOT Boost_FOUND)
  message(SEND_ERROR "Boost is required")
endif()

message(STATUS "Boost libs in ${Boost_LIBRARY_DIRS}")

find_package(PythonInterp 2.7)

if(PYTHONINTERP_FOUND)
  find_package(PythonLibs ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}.${PYTHON_VERSION_PATCH})

  execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pyconfig.py
    OUTPUT_VARIABLE _NUMPY_DIRS
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  find_path(NUMPY_INCLUDE_DIRS numpy/numpyconfig.h
    PATHS ${_NUMPY_DIRS}
    NO_DEFAULT_PATH
  )
  message(STATUS "Numpy headers ${NUMPY_INCLUDE_DIRS}")
endif()

if(PYTHONLIBS_FOUND)
  find_library(Boost_PYTHON_LIBRARY
    "boost_python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}"
    PATHS ${Boost_LIBRARY_DIRS}
    NO_DEFAULT_PATH
  )

  message(STATUS "Boost python ${Boost_PYTHON_LIBRARY}")
endif()

add_subdirectory(src)

if(USE_PYTHON AND NUMPY_INCLUDE_DIRS AND Boost_PYTHON_LIBRARY)
  message(STATUS "Building python support")
  add_subdirectory(python/uscsi)
else()
  message(STATUS "Omitting python support")
endif()

enable_testing()
add_subdirectory(test)
