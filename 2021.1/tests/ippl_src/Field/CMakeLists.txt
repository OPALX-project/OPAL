set (_SRCS
    BCond.cpp
    Cell.cpp
    DataConnect.cpp
    Eureka.cpp
    Field.cpp
    FieldDebug.cpp
    FieldDebug2.cpp
    FunctionFace.cpp
    Periodic.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
