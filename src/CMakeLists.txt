include_directories(
  ${Boost_INCLUDE_DIRS}
)

add_library(uscsi_core STATIC
  scsi/util.h
  scsi/base.h
  scsi/linear.h
  scsi/state/vector.h
  scsi/state/matrix.h
  base.cpp
  linear.cpp
)
