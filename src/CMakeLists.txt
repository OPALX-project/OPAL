set (DKS_VERSION_OPAL \"${DKS_VERSION}\")

add_definitions (
    "${IPPL_CMAKE_CXX_FLAGS} -DPARALLEL_IO ${CMAKE_CXX_FLAGS}"
    )

if (ENABLE_SAAMG_SOLVER)
    set (CMAKE_CXX_FLAGS
        "-DHAVE_SAAMG_SOLVER ${CMAKE_CXX_FLAGS}"
        )
    set (SAAMG_SOLVER_LIBS
        ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
        ${Fortran_LIBS}
        )
else ()
    set (SAAMG_SOLVER_LIBS "")
endif ()

if (ENABLE_AMR)
    set(AMR_LIBS
        ${MPI_Fortran_LIBRARIES}
        ${Fortran_LIBS}
        ${AMREX_EXTRA_Fortran_LINK_LINE}
        ${AMREX_LIBRARIES}
        
    )
else (ENABLE_AMR)
    set(AMR_LIBS "")
endif (ENABLE_AMR)

if (ENABLE_AMR AND ENABLE_AMR_MG_SOLVER)
    set (CMAKE_CXX_FLAGS
        "-DHAVE_AMR_MG_SOLVER ${CMAKE_CXX_FLAGS}"
    )
    set (AMR_MG_SOLVER_LIBS
        ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
else (ENABLE_AMR AND ENABLE_AMR_MG_SOLVER)
    set (AMR_MG_SOLVER_LIBS "")
endif(ENABLE_AMR AND ENABLE_AMR_MG_SOLVER)

configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/OPALconfig.h)
configure_file(revision.h.in ${CMAKE_CURRENT_SOURCE_DIR}/OPALrevision.h)

set (OPAL_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set (OPAL_SRCS )
macro (ADD_OPAL_SOURCES)
    file (RELATIVE_PATH _relPath "${OPAL_SRC_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND OPAL_SRCS "${_relPath}/${_src}")
        else()
            list (APPEND OPAL_SRCS "${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SRCS to parent directory
        set (OPAL_SRCS ${OPAL_SRCS} PARENT_SCOPE)
    endif()
endmacro()

macro (ADD_OPAL_HEADERS )
    file (RELATIVE_PATH _relPath "${OPAL_SRC_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_hdr ${ARGN})
        if (_relPath)
            list (APPEND OPAL_HDRS "${_relPath}/${_hdr}")
        else ()
            list (APPEND OPAL_HDRS "${_hdr}")
        endif ()
    endforeach ()
    if (_relPath)
        # propagate HDRS to parent directory
        set (OPAL_HDRS ${OPAL_HDRS} PARENT_SCOPE)
    endif ()
endmacro ()

ADD_OPAL_SOURCES(opal.cpp)
ADD_OPAL_SOURCES(changes.cpp)

add_subdirectory (AbstractObjects)
add_subdirectory (Algebra)
add_subdirectory (Algorithms)
add_subdirectory (Aperture)
add_subdirectory (Attributes)
add_subdirectory (BasicActions)
add_subdirectory (Classic)
add_subdirectory (Distribution)
add_subdirectory (Editor)
add_subdirectory (Elements)
add_subdirectory (Expressions)
add_subdirectory (Lines)
add_subdirectory (Match)
add_subdirectory (OpalConfigure)
add_subdirectory (OpalParser)
add_subdirectory (Optimize)
add_subdirectory (PhysicsActions)
add_subdirectory (Sample)
add_subdirectory (Solvers)
add_subdirectory (Steppers)
add_subdirectory (Structure)
add_subdirectory (Tables)
add_subdirectory (Track)
add_subdirectory (Utilities)
add_subdirectory (ValueDefinitions)

if (ENABLE_AMR)
    add_subdirectory (Amr)
endif ()

include_directories ( BEFORE
    ${OPAL_SRC_DIR}
    ${CLASSIC_SRC_DIR}
    ${IPPL_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${GSL_INCLUDE_DIR}
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
    ${IPPL_SOURCE_DIRS}
    ${AMREX_INCLUDE_DIR}
    ${GTEST_INCLUDE_DIRS}
    ${OPTP_INCLUDE_DIR}
    )

link_directories (
    ${IPPL_LIBRARY_DIR}
    ${OPTP_LIBRARY_DIR}
    ${Trilinos_LIBRARY_DIRS}
    ${Trilinos_TPL_LIBRARY_DIRS}
    ${AMREX_LIBRARY_DIR}
    )

if (ENABLE_DKS)
    include_directories (${DKS_INCLUDE_DIR})
    link_directories (${DKS_LIBRARY_DIR})
endif ()

#
# we cannot use the target name 'OPAL' for 'libOPAL'. On systems with
# case-insensitive filesystems (like macOS) this conflicts with the
# target 'opal' for the OPAL binary.
#
add_library (libOPAL ${OPAL_SRCS})
set_target_properties (libOPAL PROPERTIES OUTPUT_NAME OPAL)

if (ENABLE_DKS)
    target_link_libraries (libOPAL
        dks
        ${DKS_CUDA_LIBS})
endif ()

target_link_libraries (libOPAL
    ${SAAMG_SOLVER_LIBS}
    ${AMR_MG_SOLVER_LIBS}
    ${AMR_LIBS}
    ${MPI_CXX_LIBRARIES}
    )

add_executable (opal Main.cpp)
target_link_libraries( opal
    libOPAL
    ${OPTP_LIBS}
    ${OPTP_LIBRARY}
    ${IPPL_LIBRARY}
    ${GSL_LIBRARY}
    ${GSL_CBLAS_LIBRARY}
    ${H5Hut_LIBRARY}
    ${HDF5_LIBRARIES}
    ${Boost_LIBRARIES}
    m
    z
    ${CMAKE_DL_LIBS}
    )

if (ENABLE_OPAL_FEL)
    include_directories (${MITHRA_INCLUDE_DIR})
    #target_include_directories(opal PRIVATE ${MITHRA_INCLUDE_DIR})
endif()

install (TARGETS ${TEST_EXE} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
install (TARGETS opal RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
install (TARGETS libOPAL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

set (HDRS
    OPALconfig.h
    opal.h
    OPALrevision.h
    changes.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/src")

# vi: set et ts=4 sw=4 sts=4:
# Local Variables:
# mode:cmake
# cmake-tab-width: 4
# indent-tabs-mode:nil
# End:

#Turn on gcov
#set(CMAKE_CXX_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
#set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)
