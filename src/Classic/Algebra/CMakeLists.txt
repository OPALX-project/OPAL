set (_SRCS
    ComplexEigen.cpp
    DoubleEigen.cpp
    DynamicFixedPoint.cpp
    NormalForm.cpp
    StaticFixedPoint.cpp
    TpsComplex.cpp
    TpsData.cpp
    TpsDouble.cpp
    TpsMonomial.cpp
    TpsTps.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    Array1D.h
    Array2D.h
    ComplexEigen.h
    DoubleEigen.h
    DynamicFixedPoint.h
    LieMap.h
    LUMatrix.h
    Matrix.h
    NormalForm.h
    SliceIterator.h
    StaticFixedPoint.h
    TpsData.h
    Tps.h
    Tps.hpp
    TpsMath.h
    TpsMonomial.h
    TpsSubstitution.h
    Vector.h
    Vps.h
    Vps.hpp
    VpsInvMap.h
    VpsMap.h
    )

install (FILES ${HDRS}  DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algebra")
# Some source files included in header files
install (FILES ${_SRCS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algebra")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
