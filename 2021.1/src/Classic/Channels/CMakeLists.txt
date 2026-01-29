set (_SRCS
    Channel.cpp
    ConstChannel.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})


set (HDRS
    Channel.h
    ConstChannel.h
    DirectChannel.h
    IndexedChannel.h
    IndirectChannel.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Channels")
