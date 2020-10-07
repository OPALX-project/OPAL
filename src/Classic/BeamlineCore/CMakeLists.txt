set (_SRCS
    BeamStrippingRep.cpp
    CCollimatorRep.cpp
    CorrectorRep.cpp
    CyclotronRep.cpp
    DegraderRep.cpp
    DriftRep.cpp
    FlexibleCollimatorRep.cpp
    MarkerRep.cpp
    MonitorRep.cpp
    MultipoleRep.cpp
    Octupole.cpp
    ProbeRep.cpp
    Quadrupole.cpp
    RBendRep.cpp
    RFCavityRep.cpp
    SBendRep.cpp
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
    YCorrectorRep.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    BeamStrippingRep.h
    CCollimatorRep.h
    CorrectorRep.h
    CyclotronRep.h
    DegraderRep.h
    DriftRep.h
    FlexibleCollimatorRep.h
    MarkerRep.h
    MonitorRep.h
    MultipoleRep.h
    Octupole.h
    ProbeRep.h
    Quadrupole.h
    RBendRep.h
    RFCavityRep.h
    SBendRep.h
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
    YCorrectorRep.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/BeamlineCore")
