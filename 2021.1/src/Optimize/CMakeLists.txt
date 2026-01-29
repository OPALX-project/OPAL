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
