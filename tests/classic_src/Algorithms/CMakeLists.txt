set (_SRCS
    PolynomialTimeDependenceTest.cpp
    SplineTimeDependenceTest.cpp
    DistributionMomentsTest.cpp
    DistributionMomentsTestFixture.cpp
    PartBunchTest.cpp
    SinusoidalTimeDependenceTest.cpp
    TimeDependenceCacheTest.cpp
    CoreFitSCTest.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
