set (_SRCS
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    BorisPusher.h
    LF2.h
    LF2.hpp
    RK4.h
    RK4.hpp
    Stepper.h
    Steppers.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Steppers")
