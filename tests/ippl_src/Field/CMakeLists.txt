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

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End: