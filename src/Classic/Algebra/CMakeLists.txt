set (_SRCS
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
    LUMatrix.h
    Matrix.h
    SliceIterator.h
    TpsData.h
    Tps.h
    Tps.hpp
    TpsMath.h
    TpsMonomial.h
    TpsSubstitution.h
    Vector.h
    Vps.h
    Vps.hpp
    VpsMap.h
    )

install (FILES ${HDRS}  DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algebra")
# Some source files included in header files
install (FILES ${_SRCS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algebra")
