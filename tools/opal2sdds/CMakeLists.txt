cmake_minimum_required (VERSION 2.8.10)
project (OPAL2SDDS)
set (OPAL2SDDS_VERSION_MAJOR 0)
set (OPAL2SDDS_VERSION_MINOR 1)

configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

add_definitions (-DPARALLEL_IO)

include_directories (
    ${SDDS_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
)

link_directories (
    ${SDDS_LIBRARY_DIR}
)

set (SDDS_LIBS
    mdbcommon
    SDDS1
    mdblib
    lzma
    ${H5Hut_LIBRARY}
    ${HDF5_LIBRARIES}
    z
    ${MPI_LIBRARIES}
    dl
)

add_executable (opal2sdds main.cpp)
target_link_libraries (opal2sdds ${SDDS_LIBS})

INSTALL(TARGETS opal2sdds RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")

# vi: set et ts=4 sw=4 sts=4:
  
# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
