set (_SRCS
    DragtFinnMap.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    DragtFinnMap.h
    DragtFinnNormalForm.h
    FArray1D.h
    FArray2D.h
    FDoubleEigen.h
    FDynamicFP.h
    FLieGenerator.h
    FLieGenerator.hpp
    FLUMatrix.h
    FMatrix.h
    FMonomial.h
    FNormalForm.h
    FSlice.h
    FStaticFP.h
    FTpsData.h
    FTps.h
    FTps.hpp
    FTpsMath.h
    FVector.h
    FVps.h
    FVps.hpp
    LinearFun.h
    LinearFun.hpp
    LinearMap.h
    LinearMap.hpp
    LinearMath.h
    Taylor.h
    Taylor.hpp
    TransportFun.h
    TransportFun.hpp
    TransportMap.h
    TransportMap.hpp
    TransportMath.h
    )

install (FILES ${HDRS}  DESTINATION "${CMAKE_INSTALL_PREFIX}/include/FixedAlgebra")
# Some source files included in header files
install (FILES ${_SRCS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/FixedAlgebra")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
