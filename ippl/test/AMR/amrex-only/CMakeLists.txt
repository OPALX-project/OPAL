file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message (STATUS "Adding test AMR found in ${_relPath}")

add_compile_options (-O3)
add_compile_options (-stdc=c++11)
   
add_definitions (-DIPPL_AMR)
add_definitions (${AMREX_DEFINES})

include_directories (
    ${CMAKE_SOURCE_DIR}/ippl/src
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/ippl/src
    ${CMAKE_SOURCE_DIR}/src/Classic/
    ${CMAKE_SOURCE_DIR}/src/Classic/Utilities/
    ${CMAKE_SOURCE_DIR}/src/Utilities/
    ${IPPL_INCLUDE_DIR}
    ${AMREX_INCLUDE_DIR}
)
    
link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${IPPL_LIBRARY_DIR}
    ${AMREX_LIBRARY_DIR}
    ${HYPRE_LIBRARY_DIR}
)

set (LIBS
    ${IPPL_LIBRARY}
    ippl
    ${MPI_CXX_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
    ${Fortran_LIBS}
    ${AMREX_EXTRA_Fortran_LINK_LINE}
    ${AMREX_LIBRARIES}
    m
    z
)

if (${AMREX_DIM} EQUAL 3)
    add_executable (testDeposition
        testDeposition.cpp)
    
    add_executable (AssignMultiLevelDensity
        AssignMultiLevelDensity.cpp)

    #add_executable( testGeometry
    #   testGeometry.cpp)
    
    add_executable( testDepositionFail
        testDepositionFail.cpp)
    
    add_executable( testSolverFail
        testSolverFail.cpp
        ../Solver.cpp)
    
    # 4. March 2017,
    # http://stackoverflow.com/questions/13429656/cmake-how-to-copy-contents-of-a-directory-into-build-directory-after-make
    add_custom_command(TARGET testDeposition testDepositionFail POST_BUILD
        COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/inputs $<TARGET_FILE_DIR:testDeposition>)

    target_link_libraries (testDeposition ${LIBS}
        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
                        
    target_link_libraries (AssignMultiLevelDensity ${LIBS}
        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    #target_link_libraries (testGeometry ${LIBS}
    #   ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testDepositionFail ${LIBS}
        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})
    
    target_link_libraries (testSolverFail ${LIBS}
        ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

    set (EXECS
        testDeposition
        AssignMultiLevelDensity
        testDepositionFail
        testSolverFail
        )

    install (TARGETS ${EXECS} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/libexec/amr/amrex-only")
endif ( ${AMREX_DIM} EQUAL 3)
