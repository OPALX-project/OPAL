set (_SRCS
    AbsFileStream.cpp
    FileStream.cpp
    SimpleStatement.cpp
    Statement.cpp
    StringStream.cpp
    TerminalStream.cpp
    Token.cpp
    TokenStream.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    AbsFileStream.h
    FileStream.h
    Parser.h
    SimpleStatement.h
    Statement.h
    StringStream.h
    TerminalStream.h
    Token.h
    TokenStream.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Parser")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
