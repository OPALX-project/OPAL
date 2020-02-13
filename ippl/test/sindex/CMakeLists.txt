file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message (STATUS "Adding index test found in ${_relPath}")

include_directories (
    ${CMAKE_SOURCE_DIR}/ippl/src
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
)

set (IPPL_LIBS ippl)
set (COMPILE_FLAGS ${OPAL_CXX_FLAGS})

if (ENABLE_DKS)
   include_directories (${DKS_INCLUDE_DIR})
   link_directories (${DKS_LIBRARY_DIR})
   SET (IPPL_LIBS ${IPPL_LIBS} ${DKS_LIBRARY_DIR}/libdks.a)
endif ()

add_executable (TestSIndex TestSIndex.cpp)
target_link_libraries (
    TestSIndex
    ${IPPL_LIBS}
    ${MPI_CXX_LIBRARIES}
    boost_timer
)

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End: