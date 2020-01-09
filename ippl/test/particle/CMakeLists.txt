file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message (STATUS "Adding test FFT found in ${_relPath}")

include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/ippl/src
    ${CLASSIC_SOURCE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
)

add_compile_options (-std=c++11)
add_definitions (-DPARALLEL_IO)

set (IPPL_LIBS ippl)
set (COMPILE_FLAGS ${OPAL_CXX_FLAGS})

if (ENABLE_DKS)
   
   include_directories (${DKS_INCLUDE_DIR})
   link_directories (${DKS_LIBRARY_DIR})
   #SET (IPPL_LIBS ${IPPL_LIBS} dks)
   SET (IPPL_LIBS ${IPPL_LIBS} ${DKS_LIBRARY_DIR}/libdks.a)

endif ()

set_source_files_properties(SOURCE PIC3d.cpp p3m3d.cpp
    PROPERTY COMPILE_FLAGS
    -fno-tree-vrp)

#add_executable (PIC3d PIC3d.cpp)
#target_link_libraries (PIC3d ${IPPL_LIBS})

add_executable (p3m3dHeating p3m3dHeating.cpp)
target_link_libraries (
    p3m3dHeating
    ${IPPL_LIBS}
    ${H5Hut_LIBRARY}
    ${HDF5_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    boost_timer
)

#add_executable (p3m3d p3m3d.cpp)
#target_link_libraries (p3m3d ${IPPL_LIBS})

#add_executable (p3m3dleapfrog p3m3dleapfrog.cpp)
#target_link_libraries (p3m3dleapfrog ${IPPL_LIBS})

#add_executable (p3m3dTwoStream p3m3dTwoStream.cpp)
#target_link_libraries (p3m3dTwoStream ${IPPL_LIBS})

add_executable (p3m3dTwoStreamParallel p3m3dTwoStreamParallel.cpp)
target_link_libraries (
    p3m3dTwoStreamParallel
    ${IPPL_LIBS}
    ${MPI_CXX_LIBRARIES}
    boost_timer
)

add_executable (p3m3dRegressionTests p3m3dRegressionTests.cpp)
target_link_libraries (
    p3m3dRegressionTests
    ${IPPL_LIBS}
    ${MPI_CXX_LIBRARIES}
    boost_timer
)

#add_executable (p3m3dMicrobunching p3m3dMicrobunching.cpp)
#target_link_libraries (p3m3dMicrobunching ${IPPL_LIBS} ${H5Hut_LIBRARY} ${HDF5_LIBRARIES} z ${CMAKE_DL_LIBS})

#add_executable (p3m3d_periodic p3m3d_periodic.cpp)
#target_link_libraries (p3m3d_periodic ${IPPL_LIBS})

#add_executable (pm pm.cpp)
#target_link_libraries (pm ${IPPL_LIBS})

#add_executable (pp pp.cpp)
#target_link_libraries (pp ${IPPL_LIBS})

#add_executable (pp_periodic pp_periodic.cpp)
#target_link_libraries (pp_periodic ${IPPL_LIBS})

#add_executable (pm_periodic pm_periodic.cpp)
#target_link_libraries (pm_periodic ${IPPL_LIBS})

#add_executable (domaintest mini_exmpl_domain.cpp)
#target_link_libraries (domaintest ${IPPL_LIBS})

#add_executable (interpolation_test interpolation_test.cpp)
#target_link_libraries (interpolation_test ${IPPL_LIBS})

#add_executable (test-scatter-1 test-scatter-1.cpp)
#target_link_libraries (test-scatter-1 ${IPPL_LIBS})

# vi: set et ts=4 sw=4 sts=4:
  
# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
