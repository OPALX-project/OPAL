set (_SRCS
    CavityAutophaser.cpp
    Ctunes.cpp
    IndexMap.cpp
    Hamiltonian.cpp
    lomb.cpp
    NilTracker.cpp
    MapAnalyser.cpp
    MultiBunchHandler.cpp
    OrbitThreader.cpp
    ParallelCyclotronTracker.cpp
    ParallelTTracker.cpp
    StepSizeConfig.cpp
    ThickMapper.cpp
    ThickTracker.cpp
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
    lomb.h
    MapAnalyser.h
    MultiBunchHandler.h
    NilTracker.h
    OrbitThreader.h
    ParallelCyclotronTracker.h
    ParallelTTracker.h
    StepSizeConfig.h
    ThickMapper.h
    ThickTracker.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
