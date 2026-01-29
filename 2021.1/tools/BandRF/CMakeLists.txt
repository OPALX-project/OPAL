message (STATUS "Configuring BandRF")

set (BANDRF_LIBS
    ${H5Hut_LIBRARY}
    ${MPI_CXX_LIBRARIES}
)

add_executable (ascii2h5block ascii2h5block.cpp)
target_link_libraries (ascii2h5block ${BANDRF_LIBS})

install (
    TARGETS ascii2h5block
    RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
)
