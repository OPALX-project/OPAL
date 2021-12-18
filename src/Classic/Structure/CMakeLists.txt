set (_SRCS
    BoundingBox.cpp
    LossDataSink.cpp
    MeshGenerator.cpp
    PeakFinder.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    BoundingBox.h
    LossDataSink.h
    MeshGenerator.h
    PeakFinder.h
    ValueRange.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Structure")