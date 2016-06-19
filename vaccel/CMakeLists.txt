include_directories(
  ${EPICS_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../src
)

add_library(flame
  flame.cpp
  flamedev.cpp
  flameelem.cpp
  flameset.cpp
  flamemeasure.cpp
  flameorbit.cpp
)
target_link_libraries(flame
  flame_core
  ${Boost_REGEX_LIBRARY}
)

epics_install(
  PREFIX lib/epics/
  LIBS flame
  DBDS flameSupport.dbd
)

if(NEED_DEMOIOC)
  epics_add_ioc(flamedemo
    DBDS flameSupport.dbd
    LIBS flame
    INSTALL_PREFIX lib/epics/
  )
endif()
