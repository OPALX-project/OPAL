file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message (STATUS "Adding test AMR found in ${_relPath}")

add_compile_options (-O3)
add_compile_options (-std=c++11)

add_definitions (-DBL_USE_FORTRAN_MPI)
add_definitions (-DUNIQUE_PTR)
add_definitions (-DHAVE_AMR_MG_SOLVER)
add_definitions (-DUSEH5FEDV2)
add_definitions (-DPARALLEL_IO)
add_definitions (-DIPPL_AMR)
add_definitions (${AMREX_DEFINES})

add_compile_options(-Wno-deprecated-declarations)

include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/ippl/src
    ${CMAKE_SOURCE_DIR}/src/Classic/
    ${CMAKE_SOURCE_DIR}/src/Classic/Utilities/
    ${CMAKE_SOURCE_DIR}/src/Utilities/
    ${IPPL_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${AMREX_INCLUDE_DIR}
    ${HYPRE_INCLUDE_DIR}
    ${BOOST_INCLUDE_DIR}
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIR}
)
    
link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../ippl-amr/
    ${CMAKE_SOURCE_DIR}/src
    ${IPPL_LIBRARY_DIR}
    ${AMREX_LIBRARY_DIR}
    ${HYPRE_LIBRARY_DIR}
    ${BOOST_LIBRARY_DIR}
    ${Trilinos_LIBRARY_DIRS}
    ${Trilinos_TPL_LIBRARY_DIRS}
    ${GSL_LIBRARY_DIR}
)

set (TESTLIBS amropal multigrid)

set (EXECS "")

if (${AMREX_DIM} EQUAL 2)
    add_executable (testRestrictionMatrix       testRestrictionMatrix.cpp)
    add_executable (testInterpolationMatrix     testInterpolationMatrix.cpp)
    add_executable (testBoundaryMatrix          testBoundaryMatrix.cpp)
    add_executable (testSmootherMatrix          testSmootherMatrix.cpp)
    add_executable (testPoissonMatrix           testPoissonMatrix.cpp)
    add_executable (testNoFine                  testNoFine.cpp)
    add_executable (testWithFine                testWithFine.cpp)
    add_executable (testSpecialPoissonMatrix    testSpecialPoissonMatrix.cpp)
    add_executable (testLagrange                testLagrange.cpp)
endif ()

# add_executable (testAmrMultigrid
#     testAmrMultigrid.cpp)

# add_executable (testSolverComparison
#     testSolverComparison.cpp
#     ../Solver.cpp)

if (${AMREX_DIM} EQUAL 3)
    set (LIBS
        ${GSL_LIBRARY}
        ${GSL_CLBAS_LIBRARY}
        ${Trilinos_LIBRARIES}
        ${Trilinos_TPL_LIBRARIES}
        ${OTHER_CMAKE_EXE_LINKER_FLAGS}
        ${CMAKE_DL_LIBS})
        )
    add_executable( testCopy testCopy.cpp)
    target_link_libraries (testCopy
        ${TESTLIBS}
        ${LIBS}
        )
    
    add_executable( testBitPattern testBitPattern.cpp)
    target_link_libraries (testBitPattern
        ${OTHER_CMAKE_EXE_LINKER_FLAGS}
        ${CMAKE_DL_LIBS}
        )
    
    add_executable( trilinos2ippl trilinos2ippl.cpp)
    target_link_libraries (trilinos2ippl
        ${TESTLIBS}
        ${LIBS}
        )
    
    set (EXECS
         ${EXECS}
         testBitPattern
         trilinos2ippl
         testCopy
         )
endif ()

if (${AMREX_DIM} EQUAL 2)
    set (LIBS
        ${Trilinos_LIBRARIES}
        ${Trilinos_TPL_LIBRARIES}
        ${OTHER_CMAKE_EXE_LINKER_FLAGS}
        ${CMAKE_DL_LIBS}
        )
    target_link_libraries (testRestrictionMatrix        ${TESTLIBS} ${LIBS})
    target_link_libraries (testInterpolationMatrix      ${TESTLIBS} ${LIBS})
    target_link_libraries (testBoundaryMatrix           ${TESTLIBS} ${LIBS})
    target_link_libraries (testSmootherMatrix           ${TESTLIBS} ${LIBS})
    target_link_libraries (testPoissonMatrix            ${TESTLIBS} ${LIBS})
    target_link_libraries (testNoFine                   ${TESTLIBS} ${LIBS})
    target_link_libraries (testWithFine                 ${TESTLIBS} ${LIBS})
    target_link_libraries (testSpecialPoissonMatrix     ${TESTLIBS} ${LIBS})
    target_link_libraries (testLagrange                 ${TESTLIBS} ${LIBS})

    set( EXECS
         ${EXECS}
         testRestrictionMatrix
         testInterpolationMatrix
         testBoundaryMatrix
         testSmootherMatrix
         testPoissonMatrix
         testNoFine
         testWithFine
         testSpecialPoissonMatrix
         testLagrange
    )
endif ()

# target_link_libraries (testAmrMultigrid               ${TESTLIBS} ${LIBS})
# target_link_libraries (testSolverComparison           ${TESTLIBS} ${LIBS})
#
# set( EXECS
#     ${EXECS}
#     testAmrMultigrid
#     testSolverComparison
#     )

install (TARGETS ${EXECS} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/libexec/amr/trilinos")

# vi: set et ts=4 sw=4 sts=4:
  
# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
