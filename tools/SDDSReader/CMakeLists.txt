cmake_minimum_required (VERSION 2.8.10)
project (StandAloneSDDSReader)
set (StandAloneSDDSReader_VERSION_MAJOR 0)
set (StandAloneSDDSReader_VERSION_MINOR 1)

configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

add_definitions (-DPARALLEL_IO)

set (PARSERPATH "${CMAKE_SOURCE_DIR}/src/Structure")

include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/Classic
)

set (SRCS
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/array.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/associate.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/ast.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/column.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/data.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/description.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/file.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/include.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/parameter.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/version.cpp
  ${CMAKE_SOURCE_DIR}/src/Utilities/OpalException.cpp
  ${CMAKE_SOURCE_DIR}/src/Classic/Utilities/ClassicException.cpp
)

set (LIBS
    ${Boost_LIBRARIES}
    ${MPI_LIBRARIES}
    z
    dl
)

add_executable (SDDSReader main.cpp ${SRCS})
target_link_libraries (SDDSReader ${LIBS})

install (TARGETS SDDSReader RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")

# vi: set et ts=4 sw=4 sts=4:
  
# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
