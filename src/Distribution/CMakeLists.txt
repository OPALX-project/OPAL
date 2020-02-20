set (_SRCS
    Distribution.cpp
    LaserProfile.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    ClosedOrbitFinder.h
    Distribution.h
    Harmonics.h
    LaserProfile.h
    MapGenerator.h
    matrix_vector_operation.h
    rdm.h
    SigmaGenerator.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Distribution")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
