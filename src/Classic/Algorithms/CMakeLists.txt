set (_SRCS
    AbstractTimeDependence.cpp
    AbstractTracker.cpp
    CoordinateSystemTrafo.cpp
    DefaultVisitor.cpp
    Flagger.cpp
    PartBunch.cpp
    PartBins.cpp
    PartBinsCyc.cpp
    PartData.cpp
    OpalParticle.cpp
    PolynomialTimeDependence.cpp
    SplineTimeDependence.cpp
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
    AbstractTimeDependence.h
    AbstractTracker.h
    CoordinateSystemTrafo.h
    DefaultVisitor.h
    Flagger.h
    ListElem.h
    PartBinsCyc.h
    PartBins.h
    PartBunch.h
    PartBunchBase.h
    PartData.h
    OpalParticle.h
    PBunchDefs.h
    PolynomialTimeDependence.h
    Quaternion.h
    SplineTimeDependence.h
    Tracker.h
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
# require-final-newline: nil
# End:
