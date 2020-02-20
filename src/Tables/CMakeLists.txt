set (_SRCS
    AttList.cpp
    AttWriter.cpp
    CorrectionBase.cpp
    Eigen.cpp
    Envelope.cpp
    Insertion.cpp
    List.cpp
    MatrixCmd.cpp
    Micado.cpp
    Period.cpp
    RangeSelector.cpp
    Selector.cpp
    Survey.cpp
    TableTester.cpp
    ThreadAll.cpp
    ThreadBpm.cpp
    Twiss.cpp
    Twiss3.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    AttList.h
    AttWriter.h
    CorrectionBase.h
    Eigen.h
    Envelope.h
    Flatten.h
    Insertion.h
    List.h
    MatrixCmd.h
    Micado.h
    Period.h
    RangeSelector.h
    Selector.h
    Survey.h
    TableTester.h
    ThreadAll.h
    ThreadBpm.h
    Twiss3.h
    Twiss.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Tables")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
