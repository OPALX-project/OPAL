set (_SRCS
    OpalOffsetTest.cpp
    OpalPolynomialTimeDependenceTest.cpp
    OpalSplineTimeDependenceTest.cpp
    OpalVariableRFCavityTest.cpp
    OpalVariableRFCavityFringeFieldTest.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
