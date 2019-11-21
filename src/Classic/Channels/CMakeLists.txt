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

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
