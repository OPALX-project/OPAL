set (_SRCS
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    Physics.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Physics")
