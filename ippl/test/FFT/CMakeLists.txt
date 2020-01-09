file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message ("Adding test FFT found in ${_relPath}")
message ("MPI CXX compiler: ${MPI_CXX_COMPILER}")
message ("CXX compiler: ${CMAKE_CXX_COMPILER}")

include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/ippl/src
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
)

set (IPPL_LIBS ippl)

if (ENABLE_DKS)
   include_directories (${DKS_INCLUDE_DIR})
   link_directories (${DKS_LIBRARY_DIR})
   set (IPPL_LIBS ${IPPL_LIBS} ${DKS_LIBRARY_DIR}/libdks.a)
endif ()


add_executable (TestFFT TestFFT.cpp)
add_executable (TestRC TestRC.cpp)
#add_executable (TestRCGPU TestRCGPU.cpp)
#add_executable (TestRCMIC TestRCMIC.cpp)
add_executable (TestFFTCos TestFFTCos.cpp)

target_link_libraries (TestFFT ${IPPL_LIBS} ${MPI_CXX_LIBRARIES} boost_timer)
target_link_libraries (TestRC ${IPPL_LIBS} ${MPI_CXX_LIBRARIES} boost_timer)
#target_link_libraries (TestRCGPU ${IPPL_LIBS} ${MPI_CXX_LIBRARIES} boost_timer)
#target_link_libraries (TestRCMIC ${IPPL_LIBS} ${MPI_CXX_LIBRARIES} boost_timer)
target_link_libraries (TestFFTCos ${IPPL_LIBS} ${MPI_CXX_LIBRARIES} boost_timer)

# vi: set et ts=4 sw=4 sts=4:
  
# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
