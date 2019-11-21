file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message (STATUS "Adding test AMR found in ${_relPath}")
    
find_package (MPI REQUIRED)
message (STATUS "The MPI C++ compiler is: ${MPI_CXX_COMPILER}")

add_definitions (-DBL_USE_FORTRAN_MPI)
add_definitions (-DPARALLEL_IO)

add_compile_options (-O3)
add_compile_options (-std=c++11)
    
include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/ippl/src
    ${CMAKE_SOURCE_DIR}/src/Classic/
    ${IPPL_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${AMREX_INCLUDE_DIR}
    ${HYPRE_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
    
link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${IPPL_LIBRARY_DIR}
    ${AMREX_LIBRARY_DIR}
    ${HYPRE_LIBRARY_DIR}
)
    
# remove digits from hostname: edison03 -> edison
string (REGEX MATCH "[^0-9]*" HOSTNAME_BASE "${HOSTNAME}")
    
if (${HOSTNAME_BASE} MATCHES "edison" OR ${HOSTNAME_BASE} MATCHES "cori" OR
        ${HOSTNAME_BASE} MATCHES "daint")
    set (LIBS
        ${H5Hut_LIBRARY}
        ${HDF5_LIBRARIES}
        ${IPPL_LIBRARY}
        ippl
        m
        z
    )
else ()
    set (LIBS
        ${H5Hut_LIBRARY}
        ${HDF5_LIBRARIES}
        ${IPPL_LIBRARY}
        ippl
        ${MPI_LIBRARIES}
        mpi_mpifh
        m
        z
    )
endif ()
    
#     add_executable (testAmrPartBunch testAmrPartBunch.cpp AmrOpal.cpp AmrPartBunch.cpp ../H5Reader.cpp Distribution.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/Classic/Physics/Physics.cpp Solver.cpp)
#     add_executable (testSolver testSolver.cpp ../H5Reader.cpp Distribution.cpp Solver.cpp)
add_executable (iterative iterative.cpp)
#     add_executable (testError error.cpp AmrPartBunch.cpp ../H5Reader.cpp Distribution.cpp PoissonProblems.cpp Solver.cpp AmrOpal.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/Classic/Physics/Physics.cpp)
add_executable (testGaussian
    testGaussian.cpp
    ../AmrOpal.cpp AmrPartBunch.cpp
    ../H5Reader.cpp ../Distribution.cpp
    ../Solver.cpp)
    
add_executable (testUniform
    testUniform.cpp
    ../AmrOpal.cpp AmrPartBunch.cpp
    ../H5Reader.cpp ../Distribution.cpp
    ../Solver.cpp)
    
add_executable (testUnifSphere
    testUnifSphere.cpp
    ../AmrOpal.cpp AmrPartBunch.cpp
    ../Solver.cpp)# HypreABecLap.H HypreABecLap.cpp)
    
add_executable (testGridSolve
    testGridSolve.cpp
    ../Solver.cpp)
    
add_executable (testUnifSphereGrid
    testUnifSphereGrid.cpp
    ../AmrOpal.cpp AmrPartBunch.cpp
    ../Solver.cpp)
    
    
add_executable(testReal
    testReal.cpp
    ../H5Reader.cpp ../Distribution.cpp
    ../AmrOpal.cpp AmrPartBunch.cpp
    ../Solver.cpp)
    
add_executable(testMultiBunch
    testMultiBunch.cpp
    ../H5Reader.cpp ../Distribution.cpp
    ../AmrOpal.cpp AmrPartBunch.cpp
    ../Solver.cpp)
    
add_executable (testH5Read testH5Read.cpp ../H5Reader.cpp)

add_executable (testHdf5 testHdf5.cpp)
    
add_executable(toFile toFile.cpp ../H5Reader.cpp ../Distribution.cpp)
    
add_executable (testTracker
    testTracker.cpp
    ../AmrOpal.cpp AmrPartBunch.cpp
    ../H5Reader.cpp ../Distribution.cpp
    ../Solver.cpp)
    
add_executable (testOldPlasma
    testOldPlasma.cpp
    ../AmrOpal.cpp AmrPartBunch.cpp
    ../H5Reader.cpp ../Distribution.cpp
    ../Solver.cpp)

add_executable (testPerfectLoadBalancing
    testPerfectLoadBalancing.cpp
    ../AmrOpal.cpp AmrPartBunch.cpp
    ../H5Reader.cpp ../Distribution.cpp
    ../Solver.cpp)

add_executable (testOldPerformance
    testOldPerformance.cpp
    ../AmrOpal.cpp AmrPartBunch.cpp
    ../H5Reader.cpp ../Distribution.cpp
    ../Solver.cpp)
    
# Boxlib has a circular dependency between -lcboxlib and -lfboxlib
# --> resolve by: -lcboxlib -lfboxlib -lcboxlib

# target_link_libraries (testAmrPartBunch ${LIBS})
# target_link_libraries (testSolver ${LIBS})
# target_link_libraries (testError ${LIBS})

target_link_libraries (iterative)

set (LIBS
    ${LIBS}
    -lgfortran
    ${AMREX_LIBRARIES}
    -lcboxlib
    ${OTHER_CMAKE_EXE_LINKER_FLAGS}
    ${CMAKE_DL_LIBS}
    )
    
target_link_libraries (testH5Read               ${LIBS})
target_link_libraries (testHdf5                 ${LIBS})
target_link_libraries (testGaussian             ${LIBS})
target_link_libraries (testUniform              ${LIBS})
target_link_libraries (testUnifSphere           ${LIBS})
target_link_libraries (testGridSolve            ${LIBS})
target_link_libraries (testUnifSphereGrid       ${LIBS})
target_link_libraries (testReal                 ${LIBS})
target_link_libraries (testMultiBunch           ${LIBS})
target_link_libraries (toFile                   ${LIBS})
target_link_libraries (testTracker              ${LIBS})
target_link_libraries (testOldPlasma            ${LIBS} -lboost_system -lboost_filesystem)
target_link_libraries (testPerfectLoadBalancing ${LIBS})
target_link_libraries (testOldPerformance       ${LIBS})

# vi: set et ts=4 sw=4 sts=4:
  
# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
