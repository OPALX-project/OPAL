set (_SRCS
    List.cpp
    RangeSelector.cpp
    Selector.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    List.h
    RangeSelector.h
    Selector.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Tables")
