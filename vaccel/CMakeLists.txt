include_directories(
  ${EPICS_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

epics_add_ioc(flame
  DBDS flameSupport.dbd
  SRCS flame.cpp
       flamedev.cpp
       flameelem.cpp
       flameset.cpp
       flamemeasure.cpp
       flameorbit.cpp
  LIBS uscsi_core
       ${Boost_REGEX_LIBRARY}
  INSTALL_PREFIX lib/epics/
)
