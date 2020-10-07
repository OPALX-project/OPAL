set (_SRCS
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    OwnPtr.h
    Pointer.h
    RCObject.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/MemoryManagement")
