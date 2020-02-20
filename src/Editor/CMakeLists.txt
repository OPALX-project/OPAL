set (_SRCS
    Edit.cpp
    EditCmd.cpp
    EditCycle.cpp
    EditEnd.cpp
    EditFlatten.cpp
    EditInstall.cpp
    EditMove.cpp
    EditParser.cpp
    EditReflect.cpp
    EditRemove.cpp
    EditReplace.cpp
    EditSelect.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    EditCmd.h
    EditCycle.h
    EditEnd.h
    EditFlatten.h
    Edit.h
    EditInstall.h
    EditMove.h
    EditParser.h
    EditReflect.h
    EditRemove.h
    EditReplace.h
    EditSelect.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Editor")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
