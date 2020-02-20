set (_SRCS
    OpalTrimCoil.cpp
    TrimCoil.cpp
    TrimCoilFit.cpp
    TrimCoilBFit.cpp
    TrimCoilPhaseFit.cpp
    TrimCoilMirrored.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    OpalTrimCoil.h
    TrimCoil.h
    TrimCoilFit.h
    TrimCoilBFit.h
    TrimCoilPhaseFit.h
    TrimCoilMirrored.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/TrimCoils")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
