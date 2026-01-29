set (_SRCS
    LossDataSink.cpp
    MeshGenerator.cpp
    PeakFinder.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    LossDataSink.h
    MeshGenerator.h
    PeakFinder.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Structure")
