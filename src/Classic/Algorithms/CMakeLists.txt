set (_SRCS
    AbstractMapper.cpp
    AbstractTimeDependence.cpp
    AbstractTracker.cpp
    CoordinateSystemTrafo.cpp
    DefaultVisitor.cpp
    Flagger.cpp
    IdealMapper.cpp
    LinearMapper.cpp
    MapIntegrator.cpp
    Mapper.cpp
    OrbitTracker.cpp
    PartBunch.cpp
    PartBins.cpp
    PartBinsCyc.cpp
    PartData.cpp
    OpalParticle.cpp
    PolynomialTimeDependence.cpp
    SplineTimeDependence.cpp
    Surveyor.cpp
    ThinMapper.cpp
    ThinTracker.cpp
    TrackIntegrator.cpp
    Tracker.cpp
    Quaternion.cpp
    )

if (ENABLE_AMR)
    list (APPEND _SRCS AmrPartBunch.cpp)
endif ()

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    AbstractMapper.h
    AbstractTimeDependence.h
    AbstractTracker.h
    CoordinateSystemTrafo.h
    DefaultVisitor.h
    Flagger.h
    IdealMapper.h
    LinearMapper.h
    ListElem.h
    MapIntegrator.h
    Mapper.h
    OrbitTracker.h
    PartBinsCyc.h
    PartBins.h
    PartBunch.h
    PartBunchBase.h
    PartData.h
    OpalParticle.h
    PBunchDefs.h
    PolynomialTimeDependence.h
    Quaternion.h
    rbendmap.h
    SplineTimeDependence.h
    Surveyor.h
    ThinMapper.h
    ThinTracker.h
    Tracker.h
    TrackIntegrator.h
    Vektor.h
    )

if (ENABLE_AMR)
    list(APPEND HDRS AmrPartBunch.h)
endif ()

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
