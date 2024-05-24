cmake_minimum_required (VERSION 3.12)
project (MSLANG)
set (MSLANG_VERSION_MAJOR 0)
set (MSLANG_VERSION_MINOR 1)

add_definitions (-DNOCTAssert)

include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/Classic
    ${CMAKE_SOURCE_DIR}/ippl/src
    ${GSL_INCLUDE_DIR}
)

link_directories (
    ${IPPL_LIBRARY_DIR}
    ${CMAKE_BINARY_DIR}/src
    ${Boost_LIBRARY_DIRS}
)

set (MSLANG_LIBS
    libOPALstatic
    ${OPAL_LIBS}
    boost_timer
)

message (STATUS "Compiling MSLang")
add_executable (mslang mslang.cpp)
target_link_libraries (mslang ${MSLANG_LIBS})

install (TARGETS mslang RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
