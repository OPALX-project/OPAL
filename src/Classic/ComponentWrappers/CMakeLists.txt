set (_SRCS
    CorrectorWrapper.cpp
    MultipoleWrapper.cpp
    SBendWrapper.cpp
    RBendWrapper.cpp
    CyclotronWrapper.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    CorrectorWrapper.h
    CyclotronWrapper.h
    MultipoleWrapper.h
    RBendWrapper.h
    SBendWrapper.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ComponentWrappers")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
