set (_SRCS
    AbsFileStream.cpp
    FileStream.cpp
    SimpleStatement.cpp
    Statement.cpp
    StringStream.cpp
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
    Token.h
    TokenStream.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Parser")