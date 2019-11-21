set (_SRCS
    FlatWriter.cpp
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
    FlatWriter.h
    Line.h
    LineTemplate.h
    Replacer.h
    Sequence.h
    SequenceMember.h
    SequenceParser.h
    SequenceTemplate.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Lines")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
