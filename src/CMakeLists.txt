include_directories(
  ${Boost_INCLUDE_DIRS}
)

add_library(uscsi_core STATIC
  scsi/util.h
  scsi/base.h
  base.cpp

  scsi/state/vector.h
  scsi/state/matrix.h

  linear.cpp
  scsi/linear.h

  moment.cpp
  scsi/moment.h
)
