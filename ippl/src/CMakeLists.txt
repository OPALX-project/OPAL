message ("configure: make IpplInfo.h")
execute_process (
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMAND date
    OUTPUT_VARIABLE DATE_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
    )
execute_process (
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMAND uname -a
    OUTPUT_VARIABLE UNAME_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
    )
execute_process (
    COMMAND whoami
    OUTPUT_VARIABLE WHOAMI_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
    )

set (IPPL_COMPILE_ARCH \"$ENV{IPPL_ARCH}\")
set (IPPL_COMPILE_LINE \"${CMAKE_CXX_FLAGS}\")
set (IPPL_COMPILE_DATE \"${DATE_OUT}\")
set (IPPL_COMPILE_MACHINE \"${UNAME_OUT}\")
set (IPPL_COMPILE_OPTIONS \"${CMAKE_CXX_FLAGS}\")
set (IPPL_COMPILE_USER \"${WHOAMI_OUT}\")

configure_file (IpplVersions.h.in ${CMAKE_CURRENT_SOURCE_DIR}/IpplVersions.h)

if (UNUSED)
    ### print the compiler flags #############################
    set (PRETTY_PRINT_WIDTH 78)
    string (TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_STR)
    set (COMPILE_FLAGS_STR
        "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_STR}}")
    string (REGEX REPLACE " +" ";" COMPILE_FLAGS_LIST ${COMPILE_FLAGS_STR})
    
    list (LENGTH COMPILE_FLAGS_LIST NUMBER_OF_COMPILE_FLAGS)
    set (CURRENT_FLAG_INDEX 0)

    message ("")
    message ("Building IPPL with:")
    while (${CURRENT_FLAG_INDEX} LESS ${NUMBER_OF_COMPILE_FLAGS})
        set (FLAGS_PRETTY_PRINT "")
        
        list (GET COMPILE_FLAGS_LIST ${CURRENT_FLAG_INDEX} CURRENT_FLAG)
        string (LENGTH ${CURRENT_FLAG} LENGTH_INCL_NEXT)
        
        while (${LENGTH_INCL_NEXT} LESS ${PRETTY_PRINT_WIDTH})
            set (FLAGS_PRETTY_PRINT
                "${FLAGS_PRETTY_PRINT} ${CURRENT_FLAG}")
            
            math (EXPR CURRENT_FLAG_INDEX "(${CURRENT_FLAG_INDEX} + 1)")
           
            if (${CURRENT_FLAG_INDEX} EQUAL ${NUMBER_OF_COMPILE_FLAGS})
                set (LENGTH_INCL_NEXT ${PRETTY_PRINT_WIDTH})
            else (${CURRENT_FLAG_INDEX} EQUAL ${NUMBER_OF_COMPILE_FLAGS})
                list (GET COMPILE_FLAGS_LIST ${CURRENT_FLAG_INDEX} CURRENT_FLAG)
                if (CURRENT_FLAG)
                    string (LENGTH ${CURRENT_FLAG} LENGTH_CURRENT_FLAG)
                    math (EXPR LENGTH_INCL_NEXT "${LENGTH_INCL_NEXT} + 1 + ${LENGTH_CURRENT_FLAG}")
                endif ()
            endif ()
        endwhile ()
        string (STRIP ${FLAGS_PRETTY_PRINT} FLAGS_PRETTY_PRINT)
        message (STATUS "${FLAGS_PRETTY_PRINT}")
    endwhile ()
    message ("")
endif ()

set (IPPL_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
macro (add_ippl_sources)
    file (RELATIVE_PATH _relPath "${IPPL_SRC_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND IPPL_SRCS "${_relPath}/${_src}")
        else ()
            list (APPEND IPPL_SRCS "${_src}")
        endif ()
    endforeach ()
    if (_relPath)
        # propagate SRCS to parent directory
        set (IPPL_SRCS ${IPPL_SRCS} PARENT_SCOPE)
    endif ()
endmacro ()

macro (add_ippl_headers)
    file (RELATIVE_PATH _relPath "${IPPL_SRC_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_hdr ${ARGN})
        if (_relPath)
            list (APPEND IPPL_HDRS "${_relPath}/${_hdr}")
        else ()
            list (APPEND IPPL_HDRS "${_hdr}")
        endif ()
    endforeach ()
    if (_relPath)
        # propagate HDRS to parent directory
        set (IPPL_HDRS ${IPPL_HDRS} PARENT_SCOPE)
    endif ()
endmacro ()

set (IPPL_BASEDIR_HDRS
    Ippl.h
    IpplVersions.h
    )
add_ippl_headers (${IPPL_BASEDIR_HDRS})

add_subdirectory (AppTypes)
add_subdirectory (DataSource)
add_subdirectory (DomainMap)
add_subdirectory (FFT)
add_subdirectory (Field)
add_subdirectory (FieldLayout)
add_subdirectory (Index)
add_subdirectory (Meshes)
add_subdirectory (Message)
add_subdirectory (Particle)
add_subdirectory (PETE)
add_subdirectory (Region)
add_subdirectory (SubField)
add_subdirectory (SubParticle)
add_subdirectory (Utility)

if (ENABLE_AMR)
    add_subdirectory(AmrParticle)
endif ()

set (MAKE_DEF
    ../test/Makefile.def
    )

include_directories (
    BEFORE ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_library ( ippl ${IPPL_SRCS} ${IPPL_SRCS_FORT} )

if (ENABLE_DKS)
    message ( "DKS include dirs: ${DKS_INCLUDE_DIR}")
    include_directories (${DKS_INCLUDE_DIR})
    link_directories (${DKS_LIBRARY_DIR})
    target_link_libraries (
        ${TARGET_LINK_LIBRARIES} 
	ippl
	dks
        ${DKS_CUDA_LIBS}
        )
else()
    target_link_libraries(
        ippl
    )
endif ()

install (FILES ${MAKE_DEF} DESTINATION lib)
install (TARGETS ippl DESTINATION lib)
install (FILES ${IPPL_BASEDIR_HDRS} DESTINATION include)

# vi: set et ts=4 sw=4 sts=4:
# Local Variables:
# mode:cmake
# cmake-tab-width: 4
# indent-tabs-mode:nil
# End:
