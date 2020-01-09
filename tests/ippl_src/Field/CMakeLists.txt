set (_SRCS
    BCond.cpp
    Eureka.cpp
    Field.cpp
    FieldDebug.cpp
)

add_compile_options(-Wno-write-strings)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End: