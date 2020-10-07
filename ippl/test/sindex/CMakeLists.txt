file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message (STATUS "Adding index test found in ${_relPath}")

include_directories (
    ${CMAKE_SOURCE_DIR}/ippl/src
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${Boost_LIBRARY_DIRS}
)

set (IPPL_LIBS ippl)
set (COMPILE_FLAGS ${OPAL_CXX_FLAGS})

add_executable (TestSIndex TestSIndex.cpp)
target_link_libraries (
    TestSIndex
    ${IPPL_LIBS}
    ${MPI_CXX_LIBRARIES}
    boost_timer
)
