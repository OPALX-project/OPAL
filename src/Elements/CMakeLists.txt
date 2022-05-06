set (_SRCS
    OpalBend.cpp
    OpalCavity.cpp
    OpalTravelingWave.cpp
    OpalCCollimator.cpp
    OpalCyclotron.cpp
    OpalDrift.cpp
    OpalECollimator.cpp
    OpalFlexibleCollimator.cpp
    OpalDegrader.cpp
    OpalElement.cpp
    OpalHKicker.cpp
    OpalKicker.cpp
    OpalMarker.cpp
    OpalMonitor.cpp
    OpalMultipole.cpp
    OpalMultipoleT.cpp
    OpalMultipoleTStraight.cpp
    OpalMultipoleTCurvedConstRadius.cpp
    OpalMultipoleTCurvedVarRadius.cpp
    OpalOctupole.cpp
    OpalOffset/OpalLocalCylindricalOffset.cpp
    OpalOffset/OpalLocalCartesianOffset.cpp
    OpalOffset/OpalGlobalCylindricalOffset.cpp
    OpalOffset/OpalGlobalCartesianOffset.cpp
    OpalProbe.cpp
    OpalPepperPot.cpp
    OpalPolynomialTimeDependence.cpp
    OpalQuadrupole.cpp
    OpalRBend.cpp
    OpalRBend3D.cpp
    OpalRCollimator.cpp
    OpalRingDefinition.cpp
    OpalSBend.cpp
    OpalSBend3D.cpp
    OpalScalingFFAMagnet.cpp
    OpalSeptum.cpp
    OpalSextupole.cpp
    OpalSlit.cpp
    OpalSolenoid.cpp
    OpalSource.cpp
    OpalSplineTimeDependence.cpp
    OpalVacuum.cpp
    OpalVariableRFCavity.cpp
    OpalVariableRFCavityFringeField.cpp
    OpalVerticalFFAMagnet.cpp
    OpalVKicker.cpp
    OpalBeamline.cpp
    OpalStripper.cpp
    )

set (HDRS
    OpalBeamline.h
    OpalBend.h
    OpalCavity.h
    OpalCCollimator.h
    OpalCyclotron.h
    OpalDegrader.h
    OpalDrift.h
    OpalECollimator.h
    OpalElement.h
    OpalHKicker.h
    OpalKicker.h
    OpalMarker.h
    OpalMonitor.h
    OpalMultipole.h
    OpalMultipoleT.h
    OpalMultipoleTStraight.h
    OpalMultipoleTCurvedConstRadius.h
    OpalMultipoleTCurvedVarRadius.h
    OpalOctupole.h
    OpalPepperPot.h
    OpalPolynomialTimeDependence.h
    OpalProbe.h
    OpalQuadrupole.h
    OpalRBend.h
    OpalRBend3D.h
    OpalRCollimator.h
    OpalRingDefinition.h
    OpalSBend3D.h
    OpalSBend.h
    OpalScalingFFAMagnet.h
    OpalSeptum.h
    OpalSextupole.h
    OpalSlit.h
    OpalSolenoid.h
    OpalSource.h
    OpalSplineTimeDependence.h
    OpalStripper.h
    OpalTravelingWave.h
    OpalVacuum.h
    OpalVariableRFCavity.h
    OpalVariableRFCavityFringeField.h
    OpalVerticalFFAMagnet.h
    OpalVKicker.h
    OpalOffset/OpalGlobalCartesianOffset.h
    OpalOffset/OpalGlobalCylindricalOffset.h
    OpalOffset/OpalLocalCartesianOffset.h
    OpalOffset/OpalLocalCylindricalOffset.h
    )

if (ENABLE_OPAL_FEL)
    list (APPEND _SRCS OpalUndulator.cpp)
    list (APPEND HDRS OpalUndulator.h)
endif ()

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Elements/OpalOffset")
