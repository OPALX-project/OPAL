add_subdirectory(EndFieldModel)
add_subdirectory(MultipoleTFunctions)

set (_SRCS
    AlignWrapper.cpp
    AttributeSet.cpp
    BeamBeam.cpp
    BeamlineVisitor.cpp
    BeamStripping.cpp
    Bend2D.cpp
    BendBase.cpp
    CCollimator.cpp
    Component.cpp
    Corrector.cpp
    Cyclotron.cpp
    CyclotronValley.cpp
    Degrader.cpp
    Diagnostic.cpp
    Drift.cpp
    ElementBase.cpp
    ElementImage.cpp
    FlexibleCollimator.cpp
    Integrator.cpp
    Lambertson.cpp
    Marker.cpp
    Monitor.cpp
    Multipole.cpp
    MultipoleT.cpp
    MultipoleTBase.cpp
    MultipoleTStraight.cpp
    MultipoleTCurvedConstRadius.cpp
    MultipoleTCurvedVarRadius.cpp
    Offset.cpp
    ParallelPlate.cpp
    Patch.cpp
    PluginElement.cpp
    Probe.cpp
    RBend.cpp
    RBend3D.cpp
    RFCavity.cpp
    RFQuadrupole.cpp
    Ring.cpp
    SBend.cpp
    SBend3D.cpp
    ScalingFFAMagnet.cpp
    Separator.cpp
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
    AlignWrapper.h
    AttributeSet.h
    BeamBeam.h
    BeamlineVisitor.h
    BeamStripping.h
    Bend2D.h
    BendBase.h
    CCollimator.h
    Component.h
    Corrector.h
    Cyclotron.h
    CyclotronValley.h
    Degrader.h
    Diagnostic.h
    Drift.h
    ElementBase.h
    ElementImage.h
    FlexibleCollimator.h
    Integrator.h
    Lambertson.h
    Marker.h
    Monitor.h
    Multipole.h
    MultipoleT.h
    MultipoleTBase.h
    MultipoleTStraight.h
    MultipoleTCurvedConstRadius.h
    MultipoleTCurvedVarRadius.h
    Offset.h
    ParallelPlate.h
    Patch.h
    PluginElement.h
    Probe.h
    RBend.h
    RBend3D.h
    RFCavity.h
    RFQuadrupole.h
    Ring.h
    SBend3D.h
    SBend.h
    ScalingFFAMagnet.h
    SectorFieldMapComponent.h
    Separator.h
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

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
