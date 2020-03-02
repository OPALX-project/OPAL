set (_SRCS
    CavityAutophaser.cpp
    Ctunes.cpp
    IndexMap.cpp
    Hamiltonian.cpp
    LieMapper.cpp
    lomb.cpp
    NilTracker.cpp
    MapAnalyser.cpp
    MPSplitIntegrator.cpp
    MultiBunchHandler.cpp
    OrbitThreader.cpp
    ParallelCyclotronTracker.cpp
    ParallelSliceTracker.cpp
    ParallelTTracker.cpp
    StepSizeConfig.cpp
    ThickMapper.cpp
    ThickTracker.cpp
    TransportMapper.cpp
    bet/EnvelopeSlice.cpp
    bet/EnvelopeBunch.cpp
    bet/profile.cpp
    bet/math/rk.cpp
    bet/math/linfit.cpp
    bet/math/root.cpp
    bet/math/savgol.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    CavityAutophaser.h
    Ctunes.h
    Hamiltonian.h
    IndexMap.h
    LieMapper.h
    lomb.h
    MapAnalyser.h
    MPSplitIntegrator.h
    MultiBunchHandler.h
    NilTracker.h
    OrbitThreader.h
    ParallelCyclotronTracker.h
    ParallelSliceTracker.h
    ParallelTTracker.h
    StepSizeConfig.h
    ThickMapper.h
    ThickTracker.h
    TransportMapper.h
    bet/EnvelopeBunch.h
    bet/EnvelopeSlice.h
    bet/profile.h
    bet/math/linfit.h
    bet/math/rk.h
    bet/math/root.h
    bet/math/savgol.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
