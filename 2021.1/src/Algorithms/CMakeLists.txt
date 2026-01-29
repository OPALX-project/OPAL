set (_SRCS
    CavityAutophaser.cpp
    Ctunes.cpp
    IndexMap.cpp
    Hamiltonian.cpp
    lomb.cpp
    MapAnalyser.cpp
    MultiBunchHandler.cpp
    OrbitThreader.cpp
    ParallelCyclotronTracker.cpp
    ParallelTTracker.cpp
    StepSizeConfig.cpp
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
    OrbitThreader.h
    ParallelCyclotronTracker.h
    ParallelTTracker.h
    StepSizeConfig.h
    ThickTracker.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")
