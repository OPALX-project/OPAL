set (_SRCS
    BDipoleField.cpp
    BMultipoleField.cpp
    ConstBField.cpp
    ConstBzField.cpp
    ConstEField.cpp
    ConstEzField.cpp
    EDipoleField.cpp
    EMField.cpp
    NullField.cpp
    StaticElectricField.cpp
    StaticMagneticField.cpp
    Fieldmap.cpp
    FM3DH5Block.cpp
    FM3DH5Block_nonscale.cpp
    FM3DMagnetoStaticH5Block.cpp
    FM3DMagnetoStatic.cpp
    FM3DMagnetoStaticExtended.cpp
    FM3DDynamic.cpp
    FM2DElectroStatic.cpp
    FM2DMagnetoStatic.cpp
    FM2DDynamic.cpp
    FM1DDynamic.cpp
    FM1DDynamic_fast.cpp
    Astra1D_fast.cpp
    Astra1DDynamic.cpp
    Astra1DDynamic_fast.cpp
    FM1DElectroStatic.cpp
    FM1DElectroStatic_fast.cpp
    Astra1DElectroStatic.cpp
    Astra1DElectroStatic_fast.cpp
    FM1DMagnetoStatic.cpp
    FM1DMagnetoStatic_fast.cpp
    Astra1DMagnetoStatic.cpp
    Astra1DMagnetoStatic_fast.cpp
    FM1DProfile1.cpp
    FM1DProfile2.cpp
    FMDummy.cpp
    SectorMagneticFieldMap.cpp
    SectorField.cpp
    )

include_directories (
    Interpolation
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_subdirectory (Interpolation)

add_opal_sources (${_SRCS})

set (HDRS
    AcceleratingField.h
    Astra1DDynamic_fast.h
    Astra1DDynamic.h
    Astra1DElectroStatic_fast.h
    Astra1DElectroStatic.h
    Astra1D_fast.h
    Astra1DMagnetoStatic_fast.h
    Astra1DMagnetoStatic.h
    BDipoleField.h
    BMultipoleField.h
    BSingleMultipoleField.h
    ConstBField.h
    ConstBzField.h
    ConstEField.h
    ConstEzField.h
    EDipoleField.h
    EMField.h
    Fieldmap.h
    Fieldmap.hpp
    FM1DDynamic_fast.h
    FM1DDynamic.h
    FM1DElectroStatic_fast.h
    FM1DElectroStatic.h
    FM1DMagnetoStatic_fast.h
    FM1DMagnetoStatic.h
    FM1DProfile1.h
    FM1DProfile2.h
    FM2DDynamic.h
    FM2DElectroStatic.h
    FM2DMagnetoStatic.h
    FM3DDynamic.h
    FM3DH5Block.h
    FM3DH5Block_nonscale.h
    FM3DMagnetoStaticH5Block.h
    FM3DMagnetoStatic.h
    FM3DMagnetoStaticExtended.h
    FMDummy.h
    NullField.h
    OscillatingField.h
    SectorField.h
    SectorMagneticFieldMap.h
    StaticElectricField.h
    StaticMagneticField.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Fields")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
