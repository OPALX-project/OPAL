set (_SRCS
    BeamStrippingPhysics.cpp
    CollimatorPhysics.cpp
    CSRWakeFunction.cpp
    CSRIGFWakeFunction.cpp
    GreenWakeFunction.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    BeamStrippingPhysics.hh
    GreenWakeFunction.hh
    CollimatorPhysics.hh
    CSRWakeFunction.hh
    ParticleMatterInteractionHandler.hh
    WakeFunction.hh
    CSRIGFWakeFunction.hh
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Solvers")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
