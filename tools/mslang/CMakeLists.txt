cmake_minimum_required (VERSION 2.8.10)
cmake_policy (SET CMP0048 NEW)
project (MSLANG)
set (MSLANG_VERSION_MAJOR 0)
set (MSLANG_VERSION_MINOR 1)

add_definitions (-DNOCTAssert)

include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/Classic
    ${CMAKE_SOURCE_DIR}/ippl/src
    ${GSL_INCLUDE_DIR}
)

link_directories (
    ${IPPL_LIBRARY_DIR}
    ${CMAKE_SOURCE_DIR}/src
)

set (MSLANG_LIBS
    libOPAL
    ${OPTP_LIBS}
    ${OPTP_LIBRARY}
    ${IPPL_LIBRARY}
    ${GSL_LIBRARY}
    ${GSL_CBLAS_LIBRARY}
    ${H5Hut_LIBRARY}
    ${HDF5_LIBRARY}
    ${Boost_LIBRARIES}
    boost_timer
    m
    z
)

message (STATUS "Compiling MSLang")
add_executable (mslang mslang.cpp)
target_link_libraries (mslang ${MSLANG_LIBS})

install (TARGETS mslang RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")

# vi: set et ts=4 sw=4 sts=4:
  
# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
