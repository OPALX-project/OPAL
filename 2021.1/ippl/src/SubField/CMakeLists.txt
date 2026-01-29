set (_SRCS
    )

set (_HDRS
    SubBareField.hpp
    SubBareField.h
    SubFieldAssign.hpp
    SubFieldAssignDefs.h
    SubFieldAssign.h
    SubField.hpp
    SubField.h
    SubFieldIter.h
    SubFieldTraits.h
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_ippl_sources (${_SRCS})
add_ippl_headers (${_HDRS})

install (FILES ${_HDRS} DESTINATION include/SubField)
