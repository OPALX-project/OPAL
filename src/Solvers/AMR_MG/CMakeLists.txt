set (_SRCS
    AmrMultiGrid.cpp
    AmrSmoother.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

add_compile_options(-Wno-deprecated-declarations)

set (HDRS
    Amesos2BottomSolver.h
    Amesos2BottomSolver.hpp
    AmrBoundary.h
    AmrDirichletBoundary.h
    AmrInterpolater.h
    AmrLagrangeInterpolater.h
    AmrLagrangeInterpolater.hpp
    AmrMultiGridCore.h
    AmrMultiGridDefs.h
    AmrMultiGrid.h
    AmrMultiGridLevel.h
    AmrMultiGridLevel.hpp
    AmrOpenBoundary.h
    AmrPCInterpolater.h
    AmrPCInterpolater.hpp
    AmrPeriodicBoundary.h
    AmrPreconditioner.h
    AmrSmoother.h
    AmrTrilinearInterpolater.h
    AmrTrilinearInterpolater.hpp
    BelosBottomSolver.h
    BelosBottomSolver.hpp
    BottomSolver.h
    Ifpack2Preconditioner.h
    Ifpack2Preconditioner.hpp
    MueLuBottomSolver.h
    MueLuBottomSolver.hpp
    MueLuPreconditioner.h
    MueLuPreconditioner.hpp
    )


install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Solvers/AMR_MG/")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
