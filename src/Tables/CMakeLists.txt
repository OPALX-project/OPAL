set (_SRCS
    AttList.cpp
    AttWriter.cpp
    Insertion.cpp
    List.cpp
    Period.cpp
    RangeSelector.cpp
    Selector.cpp
    Survey.cpp
    Twiss.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    AttList.h
    AttWriter.h
    Flatten.h
    Insertion.h
    List.h
    Period.h
    RangeSelector.h
    Selector.h
    Survey.h
    Twiss.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Tables")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
