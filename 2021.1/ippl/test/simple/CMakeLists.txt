file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message (STATUS "Adding simple tests found in ${_relPath}")

include_directories (
    ${CMAKE_SOURCE_DIR}/ippl/src
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${Boost_LIBRARY_DIRS}
)

set (IPPL_LIBS ippl)
set (COMPILE_FLAGS ${OPAL_CXX_FLAGS})

add_executable (ScalarAdvection ScalarAdvection.cpp)
target_link_libraries (
    ScalarAdvection
    ${IPPL_LIBS}
    ${MPI_CXX_LIBRARIES}
    boost_timer
)

add_executable (fermi-dirac-field fermi-dirac-field.cpp)
target_link_libraries (
    fermi-dirac-field
    ${IPPL_LIBS}
    ${MPI_CXX_LIBRARIES}
    boost_timer
)
