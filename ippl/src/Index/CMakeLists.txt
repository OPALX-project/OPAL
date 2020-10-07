set (_SRCS
    Index.cpp
    )

set (_HDRS
    IndexedSIndex.h
    Index.h
    IndexInlines.h
    LSIndex.h
    NDIndex.h
    NDIndexInlines.h
    SIndexAssign.hpp
    SIndexAssign.h
    SIndex.hpp
    SIndex.h
    SOffset.h
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

ADD_IPPL_SOURCES (${_SRCS})
ADD_IPPL_HEADERS (${_HDRS})

install (FILES ${_HDRS} DESTINATION include/Index)
