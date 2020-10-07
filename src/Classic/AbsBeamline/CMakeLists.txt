add_subdirectory(EndFieldModel)
add_subdirectory(MultipoleTFunctions)

set (_SRCS
    AttributeSet.cpp
    BeamlineVisitor.cpp
    BeamStripping.cpp
    Bend2D.cpp
    BendBase.cpp
    CCollimator.cpp
    Component.cpp
    Corrector.cpp
    Cyclotron.cpp
    Degrader.cpp
    Drift.cpp
    ElementBase.cpp
    FlexibleCollimator.cpp
    Marker.cpp
    Monitor.cpp
    Multipole.cpp
    MultipoleT.cpp
    MultipoleTBase.cpp
    MultipoleTStraight.cpp
    MultipoleTCurvedConstRadius.cpp
    MultipoleTCurvedVarRadius.cpp
    Offset.cpp
    PluginElement.cpp
    Probe.cpp
    RBend.cpp
    RBend3D.cpp
    RFCavity.cpp
    Ring.cpp
    SBend.cpp
    SBend3D.cpp
    ScalingFFAMagnet.cpp
    Septum.cpp
    Solenoid.cpp
    Source.cpp
    Stripper.cpp
    TravelingWave.cpp
    VariableRFCavity.cpp
    VariableRFCavityFringeField.cpp
    VerticalFFAMagnet.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    AttributeSet.h
    BeamlineVisitor.h
    BeamStripping.h
    Bend2D.h
    BendBase.h
    CCollimator.h
    Component.h
    Corrector.h
    Cyclotron.h
    Degrader.h
    Drift.h
    ElementBase.h
    FlexibleCollimator.h
    Marker.h
    Monitor.h
    Multipole.h
    MultipoleT.h
    MultipoleTBase.h
    MultipoleTStraight.h
    MultipoleTCurvedConstRadius.h
    MultipoleTCurvedVarRadius.h
    Offset.h
    PluginElement.h
    Probe.h
    RBend.h
    RBend3D.h
    RFCavity.h
    Ring.h
    SBend3D.h
    SBend.h
    ScalingFFAMagnet.h
    Septum.h
    Solenoid.h
    Source.h
    SpecificElementVisitor.h
    Stripper.h
    TravelingWave.h
    VariableRFCavity.h
    VariableRFCavityFringeField.h
    VerticalFFAMagnet.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/AbsBeamline")
