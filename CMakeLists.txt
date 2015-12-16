cmake_minimum_required(VERSION 2.8.2)

project("uscsi" C CXX)

option(USE_WERROR "Treat warnings as errors")
option(USE_PYTHON "Build python support" ON)
option(NEED_PYTHON "Require python support" OFF)

message(STATUS "Option PYTHON ${USE_PYTHON}")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE DEBUG) # default to debug build
endif()
if(CMAKE_COMPILER_IS_GNUCC)
  # Turn on extra GCC warnings
  add_definitions(-Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers)
  if(USE_WERROR)
    add_definitions(-Werror)
  endif()
endif()

find_package(Boost 1.42.0 COMPONENTS system)

if(NOT Boost_FOUND)
  message(SEND_ERROR "Boost is required")
endif()

find_program(LEX_EXECUTABLE flex
  DOC "Location of lexer genereator"
)
find_program(YACC_EXECUTABLE bison
  DOC "Location of parser generator"
)
message(STATUS "FLEX ${LEX_EXECUTABLE}")
message(STATUS "BISON ${YACC_EXECUTABLE}")

if(NOT LEX_EXECUTABLE OR NOT YACC_EXECUTABLE)
  message(SEND_ERROR "Missing flex or bison")
endif()

message(STATUS "Boost libs in ${Boost_LIBRARY_DIRS}")

find_package(PythonInterp 2.6)

if(PYTHONINTERP_FOUND)
  # Introspect this interpreter
  execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pyconfig.py ${CMAKE_CURRENT_BINARY_DIR}/pyconfig.cmake
    RESULT_VARIABLE _PY_CONFIG
  )
  if(_PY_CONFIG)
    message(STATUS "Failed to fetch build information for ${PYTHON_EXECUTABLE}")
  else()
    include(${CMAKE_CURRENT_BINARY_DIR}/pyconfig.cmake)
    message(STATUS "Python version ${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}.${Python_VERSION_PATCH}")
  endif()
endif()

if(Python_FOUND)
  find_path(NUMPY_INCLUDE_DIRS numpy/numpyconfig.h
    PATHS ${Python_INCLUDE_DIRS}
    NO_DEFAULT_PATH
  )
  message(STATUS "Numpy headers ${NUMPY_INCLUDE_DIRS}")

  # see if nosetests test runner is available
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -m nose.core -V
    OUTPUT_QUIET
    ERROR_QUIET
    RESULT_VARIABLE _HAVE_NOSE
  )
  if(_HAVE_NOSE)
    if(NEED_PYTHON)
      message(SEND_ERROR "Unable to build Python support")
    else()
      message(STATUS "Can't find nosetests, not running tests")
    endif()
    set(HAVE_NOSE NO)
  else()
    set(HAVE_NOSE YES)
  endif()
endif()

enable_testing()

add_subdirectory(src)

if(USE_PYTHON AND NUMPY_INCLUDE_DIRS)
  message(STATUS "Building python support")
  add_subdirectory(python/uscsi)
elseif(NEED_PYTHON)
  message(SEND_ERROR "Unable to build Python support")
else()
  message(STATUS "Omitting python support")
endif()

add_subdirectory(tools)
