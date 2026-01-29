set (_SRCS
    Line.cpp
    LineTemplate.cpp
    Replacer.cpp
    Sequence.cpp
    SequenceMember.cpp
    SequenceParser.cpp
    SequenceTemplate.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    Line.h
    LineTemplate.h
    Replacer.h
    Sequence.h
    SequenceMember.h
    SequenceParser.h
    SequenceTemplate.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Lines")
