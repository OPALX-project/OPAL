
SET (_UTIL_SRCS
    CmdArguments.cpp
    MPIHelper.cpp

    PeakReader.cpp
    ProbeReader.cpp
    SDDSParser.cpp
    SDDSParser/array.cpp
    SDDSParser/associate.cpp
    SDDSParser/ast.cpp
    SDDSParser/column.cpp
    SDDSParser/data.cpp
    SDDSParser/description.cpp
    SDDSParser/file.cpp
    SDDSParser/include.cpp
    SDDSParser/parameter.cpp
    SDDSParser/version.cpp
)

add_sources(${_UTIL_SRCS})

# vi: set et ts=4 sw=4 sts=4:
  
# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
