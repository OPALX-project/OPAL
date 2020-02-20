set (_SRCS
    Track.cpp
    TrackCmd.cpp
    TrackEnd.cpp
    TrackParser.cpp
    TrackRun.cpp
    TrackSave.cpp
    TrackStart.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    TrackCmd.h
    TrackEnd.h
    Track.h
    TrackParser.h
    TrackRun.h
    TrackSave.h
    TrackStart.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Track")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
