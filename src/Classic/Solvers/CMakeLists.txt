set (_SRCS
    BeamStrippingPhysics.cpp
    CSRIGFWakeFunction.cpp
    CSRWakeFunction.cpp
    GreenWakeFunction.cpp
    ScatteringPhysics.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    BeamStrippingPhysics.h
    CSRIGFWakeFunction.h
    CSRWakeFunction.h
    GreenWakeFunction.h
    ParticleMatterInteractionHandler.h
    ScatteringPhysics.h
    WakeFunction.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Solvers")
