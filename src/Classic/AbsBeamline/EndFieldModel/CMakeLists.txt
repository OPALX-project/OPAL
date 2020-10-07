set (_SRCS
    EndFieldModel.cpp
    Enge.cpp
    Tanh.cpp 
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    EndFieldModel.h
    Enge.h
    Tanh.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/AbsBeamline/EndFieldModel/")
