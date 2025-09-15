# python_module function is defined in src/PyOpal/CMakeLists.txt

set (VERTICALFFAMAGNET_SRCS PyVerticalFFAMagnet.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(vertical_ffa_magnet "${VERTICALFFAMAGNET_SRCS}" "elements")

set (ELEMENT_SRCS PyOpalElement.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(opal_element "${ELEMENT_SRCS}" "elements")

set (RINGDEFINITION_SRCS PyRingDefinition.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(ring_definition "${RINGDEFINITION_SRCS}" "elements")

set (LOCALCARTESIANOFFSET_SRCS PyLocalCartesianOffset.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(local_cartesian_offset "${LOCALCARTESIANOFFSET_SRCS}" "elements")

set (GLOBALCARTESIANOFFSET_SRCS PyGlobalCartesianOffset.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(global_cartesian_offset "${GLOBALCARTESIANOFFSET_SRCS}" "elements")

set (ENGE_SRCS PyEnge.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(enge "${ENGE_SRCS}" "elements")

set (ASYMMETRICENGE_SRCS PyAsymmetricEnge.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(asymmetric_enge "${ASYMMETRICENGE_SRCS}" "elements")

set (SCALINGFFAMAGNET_SRCS PyScalingFFAMagnet.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(scaling_ffa_magnet "${SCALINGFFAMAGNET_SRCS}" "elements")

set (PROBE_SRCS PyProbe.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(probe "${PROBE_SRCS}" "elements")

set (MULTIPOLET_SRCS PyMultipoleT.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(multipolet "${MULTIPOLET_SRCS}" "elements")

set (VARIABLERFCAVITY_SRCS PyVariableRFCavity.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(variable_rf_cavity "${VARIABLERFCAVITY_SRCS}" "elements")

set (VARIABLERFCAVITYFRINGEFIELD_SRCS PyVariableRFCavityFringeField.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(variable_rf_cavity_fringe_field "${VARIABLERFCAVITYFRINGEFIELD_SRCS}" "elements")

set (POLYNOMIALTIMEDEPENDENCE_SRCS PyPolynomialTimeDependence.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(polynomial_time_dependence "${POLYNOMIALTIMEDEPENDENCE_SRCS}" "elements")

set (SINUSOIDALTIMEDEPENDENCE_SRCS PySinusoidalTimeDependence.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(sinusoidal_time_dependence "${SINUSOIDALTIMEDEPENDENCE_SRCS}" "elements")

set (SPLINETIMEDEPENDENCE_SRCS PySplineTimeDependence.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(spline_time_dependence "${SPLINETIMEDEPENDENCE_SRCS}" "elements")

set (OUTPUTPLANE_SRCS PyOutputPlane.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(output_plane "${OUTPUTPLANE_SRCS}" "elements")
