set (_SRCS
    PortableBitmapReaderTest.cpp
    PortableGraymapReaderTest.cpp
    RingSectionTest.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
