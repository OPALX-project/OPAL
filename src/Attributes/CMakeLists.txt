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
    UpperCaseString.cpp
    UpperCaseStringArray.cpp
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
    UpperCaseString.h
    UpperCaseStringArray.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Attributes")
