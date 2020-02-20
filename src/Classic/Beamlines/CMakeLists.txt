set (_SRCS
    Beamline.cpp
    BeamlineGeometry.cpp
    ElmPtr.cpp
    FlaggedBeamline.cpp
    FlaggedElmPtr.cpp
    SimpleBeamline.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    BeamlineGeometry.h
    Beamline.h
    ElmPtr.h
    FlaggedBeamline.h
    FlaggedElmPtr.h
    SimpleBeamline.h
    TBeamline.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Beamlines")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
