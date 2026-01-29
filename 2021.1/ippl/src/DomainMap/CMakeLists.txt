set (_SRCS
    )

set (_HDRS
    DomainMap.hpp
    DomainMap.h
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_ippl_sources (${_SRCS})
add_ippl_headers (${_HDRS})

install (FILES ${_HDRS} DESTINATION include/DomainMap)
