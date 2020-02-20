set (_SRCS
  Beam.cpp
  OpalWake.cpp
  ParticleMatterInteraction.cpp
  SecondaryEmissionPhysics.cpp
  BoundaryGeometry.cpp
  FieldSolver.cpp
  DataSink.cpp
  ElementPositionWriter.cpp
  H5PartWrapper.cpp
  H5PartWrapperForPC.cpp
  H5PartWrapperForPS.cpp
  H5PartWrapperForPT.cpp
  MonitorStatisticsWriter.cpp
  MultiBunchDump.cpp
  IpplInfoWrapper.cpp
  H5Writer.cpp
  MemoryWriter.cpp
  LBalWriter.cpp
  StatBaseWriter.cpp
  StatWriter.cpp
  SDDSColumn.cpp
  SDDSColumnSet.cpp
  SDDSWriter.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

set (HDRS
    Beam.h
    BoundaryGeometry.h
    DataSink.h
    ElementPositionWriter.h
    FieldSolver.h
    H5PartWrapperForPC.h
    H5PartWrapperForPS.h
    H5PartWrapperForPT.h
    H5PartWrapper.h
    IpplInfoWrapper.h
    MonitorStatisticsWriter.h
    MultiBunchDump.h
    OpalWake.h
    SecondaryEmissionPhysics.h
    ParticleMatterInteraction.h
    H5Writer.h
    MemoryWriter.h
    LBalWriter.h
    StatBaseWriter.h
    StatWriter.h
    SDDSColumn.h
    SDDSColumnSet.h
    SDDSWriter.h
)

if (ENABLE_AMR)
    list(APPEND _SRCS GridLBalWriter.cpp)
    list(APPEND HDRS GridLBalWriter.h)
endif()

if (UNIX AND NOT APPLE)
    # for Linux, BSD, Solaris, Minix
    list (APPEND _SRCS MemoryProfiler.cpp)
    list (APPEND HDRS MemoryProfiler.h)
endif()

add_opal_sources (${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Structure")

# vi: set et ts=4 sw=4 sts=4:
  
# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
