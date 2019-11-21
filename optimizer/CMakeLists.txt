#-----------------------------------------------------------------------------
# options

option (ARCH_BG "Set architecture to BlueGene" OFF)

#-----------------------------------------------------------------------------
# Common sources and libs

# add_sources macro implemented by Christof Kraus
# for all common source files (Expressions and Utils)

set (_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
macro (add_sources)
    file (RELATIVE_PATH _relPath "${_SRC_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND OPTP_SRCS "${CMAKE_SOURCE_DIR}/optimizer/${_relPath}/${_src}")
        else()
            list (APPEND OPTP_SRCS "${CMAKE_SOURCE_DIR}/optimizer/${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SRCS to parent directory
        set (OPTP_SRCS ${OPTP_SRCS} PARENT_SCOPE)
    endif()
endmacro()

set (OPTP_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)

add_subdirectory (Expression)
add_subdirectory (Util)
add_subdirectory (extlib)

#-----------------------------------------------------------------------------
# Build liboptp for linking to OPAL
#

include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_library (optp ${OPTP_SRCS} ${PILOT_OPAL_SRCS})
# Set OPTP variables everywhere (childs and parents)
set (OPTP_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}" CACHE INTERNAL "" FORCE)
set (OPTP_LIBRARY "optp" CACHE INTERNAL "" FORCE)
set (OPTP_LIBS
    pthread
    m
    ${Boost_LIBRARIES}
    wfgHypervolume
    CACHE INTERNAL "" FORCE
)

#-----------------------------------------------------------------------------
# Testing

if (BUILD_OPAL_UNIT_TESTS)
    enable_testing()
    add_subdirectory(Tests)
endif ()

#-----------------------------------------------------------------------------
# Doxygen
if (ENABLE_DOXYDOC)
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        configure_file (${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
            ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY
	    )
        add_custom_target (optpilot-doc
            ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
    endif ()
endif ()

# vi: set et ts=4 sw=4 sts=4:
  
# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
