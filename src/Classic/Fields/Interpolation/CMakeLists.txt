set (_SRCS
  Interpolator3dGridTo1d.cpp
  Interpolator3dGridTo3d.cpp
  Mesh.cpp
  MMatrix.cpp
  MVector.cpp
  NDGrid.cpp
  PolynomialPatch.cpp
  PPSolveFactory.cpp
  SolveFactory.cpp
  SquarePolynomialVector.cpp
  ThreeDGrid.cpp
  TriLinearInterpolator.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    Interpolator3dGridTo1d.h
    Interpolator3dGridTo3d.h
    Mesh.h
    MMatrix.h
    MVector.h
    NDGrid.cpp
    PolynomialCoefficient.h
    PolynomialPatch.h
    PPSolveFactory.h
    SolveFactory.h
    SquarePolynomialVector.h
    ThreeDGrid.h
    TriLinearInterpolator.h
    VectorMap.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Interpolation")
