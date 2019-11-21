set (_SRCS
    AbstractFun.cpp
    AbstractVar.cpp
    ConcreteFun.cpp
    ConcreteVar.cpp
    ConstraintCmd.cpp
    LMDif.cpp
    Match.cpp
    MatchCmd.cpp
    MatchEnd.cpp
    MatchOption.cpp
    MatchParser.cpp
    Migrad.cpp
    Simplex.cpp
    VaryCmd.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    AbstractFun.h
    AbstractVar.h
    ConcreteFun.h
    ConcreteVar.h
    ConstraintCmd.h
    LMDif.h
    MatchCmd.h
    MatchEnd.h
    Match.h
    MatchLimits.h
    MatchOption.h
    MatchParser.h
    MatchState.h
    Migrad.h
    Simplex.h
    VaryCmd.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Match")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
