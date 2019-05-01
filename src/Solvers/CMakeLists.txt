set (_SRCS
    FFTBoxPoissonSolver.cpp
    FFTPoissonSolver.cpp
    P3MPoissonSolver.cpp
 )

set (HDRS
    FFTBoxPoissonSolver.h
    FFTPoissonSolver.h
    IrregularDomain.h
    P3MPoissonSolver.h
    PoissonSolver.h
)

if ( HAVE_SAAMG_SOLVER )
    list(APPEND _SRCS
         ArbitraryDomain.cpp
         BoxCornerDomain.cpp
         EllipticDomain.cpp
         MGPoissonSolver.cpp
         RectangularDomain.cpp
         )

    list(APPEND HDRS
         ArbitraryDomain.h
         BoxCornerDomain.h
         EllipticDomain.h
         MGPoissonSolver.h
         RectangularDomain.h
         )

endif (HAVE_SAAMG_SOLVER )

set(AMR_MG_DIR "")

if ( ENABLE_AMR )
    list(APPEND _SRCS BoxLibSolvers/FMGPoissonSolver.cpp AMReXSolvers/MLPoissonSolver.cpp)
    list(APPEND HDRS AmrPoissonSolver.h BoxLibSolvers/FMGPoissonSolver.h AMReXSolvers/MLPoissonSolver.h)

    if ( ENABLE_AMR_MG_SOLVER )
        add_subdirectory (AMR_MG)
        set(AMR_MG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/AMR_MG")
    endif ( ENABLE_AMR_MG_SOLVER )

endif ( ENABLE_AMR )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${AMR_MG_DIR}
)

add_opal_sources(${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Solvers")