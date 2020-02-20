set (_SRCS
    SavitzkyGolay.cpp
    Stencil.cpp
    FixedFFTLowPass.cpp
    RelativeFFTLowPass.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    Filter.h
    Filters.h
    FixedFFTLowPass.h
    RelativeFFTLowPass.h
    SavitzkyGolay.h
    Stencil.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Filters")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
