set (_SRCS
  BeamBeamRep.cpp
  BeamStrippingRep.cpp
  CCollimatorRep.cpp
  CorrectorRep.cpp
  CyclotronRep.cpp
  CyclotronValleyRep.cpp
  DegraderRep.cpp
  DriftRep.cpp
  UndulatorRep.cpp
  FlexibleCollimatorRep.cpp
  MarkerRep.cpp
  MonitorRep.cpp
  MultipoleRep.cpp
  Octupole.cpp
  ParallelPlateRep.cpp
  PatchRep.cpp
  ProbeRep.cpp
  Quadrupole.cpp
  RBendRep.cpp
  RFCavityRep.cpp
  SBendRep.cpp
  SeparatorRep.cpp
  SeptumRep.cpp
  Sextupole.cpp
  SkewOctupole.cpp
  SkewQuadrupole.cpp
  SkewSextupole.cpp
  SolenoidRep.cpp
  SourceRep.cpp
  StripperRep.cpp
  TravelingWaveRep.cpp
  XCorrectorRep.cpp
  XMonitorRep.cpp
  YCorrectorRep.cpp
  YMonitorRep.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    BeamBeamRep.h
    BeamStrippingRep.h
    CCollimatorRep.h
    CorrectorRep.h
    CyclotronRep.h
    CyclotronValleyRep.h
    DegraderRep.h
    DriftRep.h
    UndulatorRep.h	
    FlexibleCollimatorRep.h
    MarkerRep.h
    MonitorRep.h
    MultipoleRep.h
    Octupole.h
    ParallelPlateRep.h
    PatchRep.h
    ProbeRep.h
    Quadrupole.h
    RBendRep.h
    RFCavityRep.h
    SBendRep.h
    SeparatorRep.h
    SeptumRep.h
    Sextupole.h
    SingleMultipole.h
    SkewOctupole.h
    SkewQuadrupole.h
    SkewSextupole.h
    SolenoidRep.h
    SourceRep.h
    StripperRep.h
    TravelingWaveRep.h
    XCorrectorRep.h
    XMonitorRep.h
    YCorrectorRep.h
    YMonitorRep.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/BeamlineCore")