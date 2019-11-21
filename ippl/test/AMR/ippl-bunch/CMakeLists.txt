file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message (STATUS "Adding test AMR found in ${_relPath}")

add_compile_options (-O3)
add_compile_options (-std=c++11)

add_definitions (-DPARALLEL_IO)
add_definitions (-DIPPL_AMR)
add_definitions(${AMREX_DEFINES})
    
include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/ippl/src
    ${CMAKE_SOURCE_DIR}/src/Classic/
    ${CMAKE_SOURCE_DIR}/src/Utilities/
    ${IPPL_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${AMREX_INCLUDE_DIR}
    ${HYPRE_INCLUDE_DIR}
)
    
link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${IPPL_LIBRARY_DIR}
    ${AMREX_LIBRARY_DIR}
    ${HYPRE_LIBRARY_DIR}
)
    
SET (LIBS
    ${H5Hut_LIBRARY}
    ${HDF5_LIBRARIES}
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

# add_executable(testCRTP testCRTP.cpp)

add_executable( testAdapterPattern
                testAdapterPattern.cpp
                ${CMAKE_SOURCE_DIR}/src/Classic/Utilities/ClassicException.cpp
                ${CMAKE_SOURCE_DIR}/src/Utilities/OpalException.cpp
               )

# target_link_libraries (testCRTP ${LIBS}
#                         -lgfortran ${AMREX_LIBRARIES}
#                         ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

target_link_libraries (testAdapterPattern ${LIBS}
                       ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

set( EXECS
     testAdapterPattern
)

install (TARGETS ${EXECS} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/libexec/amr/ippl-bunch")

# vi: set et ts=4 sw=4 sts=4:
  
# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
