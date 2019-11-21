set (_SRCS
    Euclid3D.cpp
    Geometry.cpp
    Matrix3D.cpp
    Euclid3DGeometry.cpp
    NullGeometry.cpp
    OffsetGeometry.cpp
    PlanarArcGeometry.cpp
    VarRadiusGeometry.cpp
    RBendGeometry.cpp
    Rotation3D.cpp
    SRotatedGeometry.cpp
    StraightGeometry.cpp
    Vector3D.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    Euclid3DGeometry.h
    Euclid3D.h
    Geometry.h
    Matrix3D.h
    NullGeometry.h
    OffsetGeometry.h
    PlanarArcGeometry.h
    VarRadiusGeometry.h
    RBendGeometry.h
    Rotation3D.h
    SRotatedGeometry.h
    StraightGeometry.h
    Vector3D.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/BeamlineGeometry")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
