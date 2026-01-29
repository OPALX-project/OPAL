set (_SRCS
    NDGridTest.cpp
    PolynomialPatchTest.cpp
    PPSolveFactoryTest.cpp
    SolveFactoryTest.cpp
    SquarePolynomialVectorTest.cpp
    ThreeDGridTest.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
