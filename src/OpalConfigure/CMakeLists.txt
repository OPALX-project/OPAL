set (_SRCS
    Configure.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    Configure.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpalConfigure")
