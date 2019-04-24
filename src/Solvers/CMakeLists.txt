set (_SRCS
  P3MPoissonSolver.cpp
  FFTPoissonSolver.cpp
  FFTBoxPoissonSolver.cpp
  MGPoissonSolver.cpp
  RectangularDomain.cpp
  EllipticDomain.cpp
  ArbitraryDomain.cpp
  BoxCornerDomain.cpp
  )

set(AMR_MG_DIR "")

if ( ENABLE_AMR )
    list(APPEND _SRCS BoxLibSolvers/FMGPoissonSolver.cpp AMReXSolvers/MLPoissonSolver.cpp)
    
    if ( ENABLE_AMR_MG_SOLVER )
        add_subdirectory (AMR_MG)
        set(AMR_MG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/AMR_MG")
    endif ( ENABLE_AMR_MG_SOLVER )

endif (ENABLE_AMR )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${AMR_MG_DIR}
)

add_opal_sources(${_SRCS})

set (HDRS
    ArbitraryDomain.h
    BoxCornerDomain.h
    EllipticDomain.h
    FFTBoxPoissonSolver.h
    FFTPoissonSolver.h
    IrregularDomain.h
    MGPoissonSolver.h
    P3MPoissonSolver.h
    PoissonSolver.h
    RectangularDomain.h
)

if ( ENABLE_AMR )
    list(APPEND HDRS AmrPoissonSolver.h BoxLibSolvers/FMGPoissonSolver.h AMReXSolvers/MLPoissonSolver.h)
endif (ENABLE_AMR )


install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Solvers")