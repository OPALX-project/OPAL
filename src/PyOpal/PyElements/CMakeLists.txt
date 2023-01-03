# python_module function is defined in src/PyOpal/CMakeLists.txt

set (VERTICALFFAMAGNET_SRCS PyVerticalFFAMagnet.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(vertical_ffa_magnet "${VERTICALFFAMAGNET_SRCS}" "elements")

set (ELEMENT_SRCS PyOpalElement.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(opal_element "${ELEMENT_SRCS}" "elements")

set (RINGDEFINITION_SRCS PyRingDefinition.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(ring_definition "${RINGDEFINITION_SRCS}" "elements")

set (LOCALCARTESIANOFFSET_SRCS PyLocalCartesianOffset.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(local_cartesian_offset "${LOCALCARTESIANOFFSET_SRCS}" "elements")

set (ENGE_SRCS PyEnge.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(enge "${ENGE_SRCS}" "elements")

set (ASYMMETRICENGE_SRCS PyAsymmetricEnge.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(asymmetric_enge "${ASYMMETRICENGE_SRCS}" "elements")

set (SCALINGFFAMAGNET_SRCS PyScalingFFAMagnet.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(scaling_ffa_magnet "${SCALINGFFAMAGNET_SRCS}" "elements")

set (PROBE_SRCS PyProbe.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(probe "${PROBE_SRCS}" "elements")
