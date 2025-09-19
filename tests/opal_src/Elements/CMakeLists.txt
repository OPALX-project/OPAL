set (_SRCS
    OpalAsymmetricEngeTest.cpp
    OpalEngeTest.cpp
    OpalOffsetTest.cpp
    OpalOutputPlaneTest.cpp
    OpalPolynomialTimeDependenceTest.cpp
    OpalScalingFFAMagnetTest.cpp
    OpalSplineTimeDependenceTest.cpp
    OpalVariableRFCavityTest.cpp
    OpalVariableRFCavityFringeFieldTest.cpp
    OpalSinusoidalTimeDependenceTest.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
