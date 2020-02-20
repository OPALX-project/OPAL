set (_SRCS
    OptimizeCmd.cpp
    DVar.cpp
    Objective.cpp
    Constraint.cpp
    OpalSimulation.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    OptimizeCmd.h
    DVar.h
    Objective.h
    Constraint.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Optimize")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
