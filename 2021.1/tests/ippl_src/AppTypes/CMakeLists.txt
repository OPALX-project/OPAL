set (_SRCS
    Chsr-2.cpp
    Operators.cpp
    Tensor.cpp
    Vektor.cpp
    VektorTz.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
