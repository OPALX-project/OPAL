project (IPPL CXX)
set (IPPL_VERSION_MAJOR 1)
set (IPPL_VERSION_MINOR 1.4)

set (IPPL_VERSION_NAME "V${IPPL_VERSION_MAJOR}.${IPPL_VERSION_MINOR}")

set (IPPL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    add_definitions (-DNOPAssert)
endif ()

set (IPPL_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src" PARENT_SCOPE)
set (IPPL_LIBRARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src" PARENT_SCOPE)
set (IPPL_LIBRARY "ippl" PARENT_SCOPE)

add_subdirectory (src)

option (ENABLE_IPPLTESTS "IPPL Tests" OFF)
if (ENABLE_IPPLTESTS)
    add_subdirectory (test)
endif ()

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake )

install (
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake
    DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/${PROJECT_NAME}"
    RENAME ${PROJECT_NAME}Config.cmake
    )

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:

