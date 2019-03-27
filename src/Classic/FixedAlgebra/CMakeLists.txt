set (_SRCS
  DragtFinnMap.cpp
  FLieGenerator.cpp
  FTps.cpp
  FVps.cpp
  LinearFun.cpp
  LinearMap.cpp
  Taylor.cpp
  TransportFun.cpp
  TransportMap.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    DragtFinnMap.h
    DragtFinnNormalForm.h
    FArray1D.h
    FArray2D.h
    FComplexEigen.h
    FDoubleEigen.h
    FDynamicFP.h
    FLieGenerator.h
    FLUMatrix.h
    FMatrix.h
    FMonomial.h
    FNormalForm.h
    FSlice.h
    FStaticFP.h
    FTpsData.h
    FTps.h
    FTpsMath.h
    FVector.h
    FVps.h
    LinearFun.h
    LinearMap.h
    LinearMath.h
    Taylor.h
    TransportFun.h
    TransportMap.h
    TransportMath.h
)

install (FILES ${HDRS}  DESTINATION "${CMAKE_INSTALL_PREFIX}/include/FixedAlgebra")
# Some source files included in header files
install (FILES ${_SRCS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/FixedAlgebra")
