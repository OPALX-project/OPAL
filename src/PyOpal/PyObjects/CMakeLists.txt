set (PARSER_SRCS PyParser.cpp $<TARGET_OBJECTS:libPyOPALPyCoreNoGlobals>)
python_module(parser "${PARSER_SRCS}" objects)

set (FIELD_SRCS PyField.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(field "${FIELD_SRCS}" objects)

set (LINE_SRCS PyLine.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(line "${LINE_SRCS}" objects)

set (TRACKRUN_SRCS PyTrackRun.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(track_run "${TRACKRUN_SRCS}" objects)

set (BEAM_SRCS PyBeam.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(beam "${BEAM_SRCS}" objects)

set (DIST_SRCS PyDistribution.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(distribution "${DIST_SRCS}" objects)

set (FIELDSOLVER_SRCS PyFieldSolver.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(field_solver "${FIELDSOLVER_SRCS}" objects)

set (TRACKCMD_SRCS PyTrackCmd.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(track "${TRACKCMD_SRCS}" objects)

set (OPTION_SRCS PyOption.cpp $<TARGET_OBJECTS:libPyOPALPyCore>)
python_module(option "${OPTION_SRCS}" objects)
