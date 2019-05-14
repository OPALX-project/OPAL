set (_SRCS
  Beam.cpp
  OpalWake.cpp
  ParticleMatterInteraction.cpp
  SecondaryEmissionPhysics.cpp
  BoundaryGeometry.cpp
  FieldSolver.cpp
  DataSink.cpp
  H5PartWrapper.cpp
  H5PartWrapperForPC.cpp
  H5PartWrapperForPS.cpp
  H5PartWrapperForPT.cpp
#  SDDSParser.cpp
  MultiBunchDump.cpp
  OpalInputInterpreter.cpp
  IpplInfoWrapper.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

set (HDRS
    Beam.h
    BoundaryGeometry.h
    DataSink.h
    FieldSolver.h
    H5PartWrapperForPC.h
    H5PartWrapperForPS.h
    H5PartWrapperForPT.h
    H5PartWrapper.h
    IpplInfoWrapper.h
    MultiBunchDump.h
    OpalInputInterpreter.h
    OpalWake.h
    SecondaryEmissionPhysics.h
    ParticleMatterInteraction.h
)

if(UNIX AND NOT APPLE)
    # for Linux, BSD, Solaris, Minix
    list(APPEND _SRCS MemoryProfiler.cpp)
    list(APPEND HDRS MemoryProfiler.h)
endif()

add_opal_sources(${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Structure")
