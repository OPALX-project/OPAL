set (_SRCS
    Attributes.cpp
    Bool.cpp
    BoolArray.cpp
    Place.cpp
    Range.cpp
    Real.cpp
    RealArray.cpp
    Reference.cpp
    opalstr.cpp
    StringArray.cpp
    TableRow.cpp
    TokenList.cpp
    TokenListArray.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    Attributes.h
    BoolArray.h
    Bool.h
    opalstr.h
    Place.h
    Range.h
    RealArray.h
    Real.h
    Reference.h
    StringArray.h
    TableRow.h
    TokenListArray.h
    TokenList.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Attributes")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
