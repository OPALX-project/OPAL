set (_SRCS
    Euclid3D.cpp
    Geometry.cpp
    Matrix3D.cpp
    Euclid3DGeometry.cpp
    NullGeometry.cpp
    PlanarArcGeometry.cpp
    VarRadiusGeometry.cpp
    RBendGeometry.cpp
    Rotation3D.cpp
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
    PlanarArcGeometry.h
    VarRadiusGeometry.h
    RBendGeometry.h
    Rotation3D.h
    StraightGeometry.h
    Vector3D.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/BeamlineGeometry")
