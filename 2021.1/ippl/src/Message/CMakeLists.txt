set (_SRCS
    CRC.cpp
    CommCreator.cpp
    CommMPI.cpp
    Communicate.cpp
    Format.cpp
    MessageFunctions.cpp
    MsgBuffer.cpp
    )

set (_HDRS
    CRC.h
    CommCreator.h
    CommMPI.h
    Communicate.h
    DataTypes.h
    Format.h
    GlobalComm.h
    GlobalComm.hpp
    Message.h
    Message.hpp
    MsgBuffer.h
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