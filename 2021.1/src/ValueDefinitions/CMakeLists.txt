set (_SRCS
    BoolConstant.cpp
    RealConstant.cpp
    RealVariable.cpp
    RealVector.cpp
    StringConstant.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    BoolConstant.h
    RealConstant.h
    RealVariable.h
    RealVector.h
    StringConstant.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ValueDefinitions")
