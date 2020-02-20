set (_SRCS
    Call.cpp
    Dump.cpp
    DumpFields.cpp
    DumpEMFields.cpp
    Echo.cpp
    Help.cpp
    Option.cpp
    Save.cpp
    Select.cpp
    Show.cpp
    Stop.cpp
    Quit.cpp
    PSystem.cpp
    System.cpp
    Title.cpp
    Value.cpp
    What.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    Call.h
    DumpFields.h
    DumpEMFields.h
    Dump.h
    Echo.h
    Help.h
    Option.h
    PSystem.h
    Quit.h
    Save.h
    Select.h
    Show.h
    Stop.h
    System.h
    Title.h
    Value.h
    What.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/BasicActions")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
