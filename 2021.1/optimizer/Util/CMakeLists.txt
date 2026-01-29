
SET (_UTIL_SRCS
    CmdArguments.cpp
    MPIHelper.cpp

    PeakReader.cpp
    ProbeReader.cpp
    ProbeHistReader.cpp
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
