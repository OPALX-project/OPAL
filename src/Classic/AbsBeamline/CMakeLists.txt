add_subdirectory(EndFieldModel)
add_subdirectory(MultipoleTFunctions)

set (_SRCS
  AlignWrapper.cpp
  AttributeSet.cpp
  BeamBeam.cpp
  BeamlineVisitor.cpp
  Bend.cpp
  BendBase.cpp
  CCollimator.cpp
  Component.cpp
  Corrector.cpp
  Cyclotron.cpp
  CyclotronValley.cpp
  Degrader.cpp
  Diagnostic.cpp
  Drift.cpp
  Undulator.cpp
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
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    AlignWrapper.h
    AttributeSet.h
    BeamBeam.h
    BeamlineVisitor.h
    Bend.h
    BendBase.h
    CCollimator.h
    Component.h
    Corrector.h
    Cyclotron.h
    CyclotronValley.h
    Degrader.h
    Diagnostic.h
    Drift.h
    Undulator.h
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
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/AbsBeamline")

#Turn on gcov
#set(CMAKE_CXX_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
#set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)
