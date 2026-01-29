project (IPPL CXX)
set (IPPL_VERSION_MAJOR 1)
set (IPPL_VERSION_MINOR 1.4)

set (IPPL_VERSION_NAME "V${IPPL_VERSION_MAJOR}.${IPPL_VERSION_MINOR}")

set (IPPL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

if (${BUILD_FOR_PRODUCTION})
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