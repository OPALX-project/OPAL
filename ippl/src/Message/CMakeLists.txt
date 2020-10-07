set (_SRCS
    CRC.cpp
    Communicate.cpp
    CommCreator.cpp
    CommMPI.cpp
    MessageFunctions.cpp
    Formatter.cpp
    )

set (_HDRS
    CommCreator.h
    CommMPI.h
    Communicate.h
    CRC.h
    DataTypes.h
    Formatter.h
    GlobalComm.hpp
    GlobalComm.h
    Message.hpp
    Message.h
    Operations.h
    TagMaker.h
    Tags.h
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_ippl_sources (${_SRCS})
add_ippl_headers (${_HDRS})

install (FILES ${_HDRS} DESTINATION include/Message)
