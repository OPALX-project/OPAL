set (_SRCS
    DipoleFieldTest.cpp
    MultipoleTTest.cpp
    OffsetTest.cpp
    PolynomialTest.cpp
    RingTest.cpp
    SBend3DTest.cpp
    ScalingFFAMagnetTest.cpp
    TrimCoilTest.cpp
    VariableRFCavityTest.cpp
    VariableRFCavityFringeFieldTest.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
