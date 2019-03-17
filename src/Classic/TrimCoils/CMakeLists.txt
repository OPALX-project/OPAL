set (_SRCS
    OpalTrimCoil.cpp
    TrimCoil.cpp
    TrimCoilFit.cpp
    TrimCoilBFit.cpp
    TrimCoilPhaseFit.cpp
    TrimCoilMirrored.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    OpalTrimCoil.h
    TrimCoil.h
    TrimCoilFit.h
    TrimCoilBFit.h
    TrimCoilPhaseFit.h
    TrimCoilMirrored.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/TrimCoils")
