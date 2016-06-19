include(FindEPICSModule)

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
  message(STATUS "Building flamedemo IOC")

  find_epics_module(
    NAME autosave
    IDFILES dbd/asSupport.dbd
    DBDS asSupport.dbd
    LIBS autosave
  )
  if(autosave_FOUND)
    message(STATUS "flamedemo will include autosave")
    list(APPEND EXTRA_DBDS ${autosave_DBDS})
    list(APPEND EXTRA_LIBS ${autosave_LIBRARIES})
  endif()

  find_epics_module(
    NAME iocstats
    IDFILES dbd/devIocStats.dbd
    DBDS devIocStats.dbd
    LIBS devIocStats
  )
  if(iocstats_FOUND)
    message(STATUS "flamedemo will include iocstats")
    list(APPEND EXTRA_DBDS ${iocstats_DBDS})
    list(APPEND EXTRA_LIBS ${iocstats_LIBRARIES})
  endif()

  find_epics_module(
    NAME caputlog
    IDFILES dbd/caPutLog.dbd
    DBDS caPutLog.dbd
    LIBS caPutLog
  )
  if(caputlog_FOUND)
    message(STATUS "flamedemo will include caputlog")
    list(APPEND EXTRA_DBDS ${caputlog_DBDS})
    list(APPEND EXTRA_LIBS ${caputlog_LIBRARIES})
  endif()


  epics_add_ioc(flamedemo
    DBDS flameSupport.dbd
         ${EXTRA_DBDS}
    LIBS flame
         ${EXTRA_LIBS}
    INSTALL_PREFIX lib/epics/
  )
endif()
