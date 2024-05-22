cmake_minimum_required (VERSION 3.12)
cmake_policy(VERSION 3.12)

project (OPAL VERSION 2022.2.99)
set (PROJECT_BUGREPORT opal@lists.psi.ch)

string (REGEX REPLACE "\\.([0-9])\\." ".0\\1." PACKAGE_VERSION ${PROJECT_VERSION})
string (REGEX REPLACE "\\.([0-9])$" ".0\\1" PACKAGE_VERSION ${PACKAGE_VERSION})
string (REGEX REPLACE "\\." "" PACKAGE_VERSION ${PACKAGE_VERSION})

message (STATUS PACKAGE_VERSION=${PACKAGE_VERSION})

if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
        )
endif ()

cmake_host_system_information (RESULT HOSTNAME QUERY HOSTNAME)

message (STATUS "Build type is: " ${CMAKE_BUILD_TYPE})
message (STATUS "Host OS System: ${CMAKE_HOST_SYSTEM}")
message (STATUS "Hostname: ${HOSTNAME}")

# Select flags.
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")
set (CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

string(TOUPPER "${CMAKE_BUILD_TYPE}" uppercase_CMAKE_BUILD_TYPE)
if (uppercase_CMAKE_BUILD_TYPE MATCHES "^(DEBUG|RELWITHDEBINFO)$")
    option (BUILD_FOR_PRODUCTION "Disable all assertion checks; Only use for production" OFF)
else()
    option (BUILD_FOR_PRODUCTION "Disable all assertion checks; Only use for production" ON)
endif()

# Resolve all library dependencies
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

find_package (MPI REQUIRED)
message (STATUS "The C++ compiler identification is: ${CMAKE_CXX_COMPILER_ID}")
message (STATUS "The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
message (STATUS "The MPI C++ compiler is: ${MPI_CXX_COMPILER}")
message (STATUS "The underlying C++ compiler is: ${CMAKE_CXX_COMPILER}")

option (ENABLE_OpenMP "Use hybrid parallelism MPI-OpenMP" OFF)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

# Disable compile time assert (used by IPPL)
add_definitions (-DNOCTAssert)

add_compile_options (-Wall)
add_compile_options (-Wunused)
add_compile_options (-Wextra)
add_compile_options (-Werror)

add_compile_options (-funroll-loops)
add_compile_options (-fstrict-aliasing)

option (WILL_BUILD_STATIC_LIBRARY "Building static library" ON)
option (BUILD_OPAL_PYTHON "Building Python interface for OPAL" OFF)
if (BUILD_OPAL_PYTHON)
    option (WILL_BUILD_SHARED_LIBRARY "Building python hooks - requires dynamic libs" ON)
    find_package (Python3 COMPONENTS Development REQUIRED)
endif ()

if (${WILL_BUILD_SHARED_LIBRARY})
    add_compile_options (-fPIE -fPIC)
    add_link_options (-pie)
    # force OPAL to link to shared libraries
    if (UNIX)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
    elseif (APPLE)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib")
    else()
        MESSAGE(FATAL_ERROR "Operating system not APPLE or UNIX - quitting")
    endif()
endif ()

# compiler dependent flags
message (STATUS "The underlying C++ compiler ID is: ${CMAKE_CXX_COMPILER_ID}")
if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel" )
    add_compile_options ("-diag-disable 383")
    add_compile_options ("-diag-disable 981")
    if (ENABLE_OpenMP)
        add_compile_options (-openmp)
        add_link_options (-openmp)
    endif ()

elseif (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    add_compile_options (-Warray-bounds-pointer-arithmetic)
    # Increase template depth for boost
    add_compile_options (-ftemplate-depth=1024)
    # get rid of linking warning with boost. Looks like that boost
    # is using this visibility setting.
    add_compile_options (-fvisibility=hidden)
    add_compile_options (-fvisibility-inlines-hidden)
    add_compile_options (-Wno-deprecated-declarations)
    add_compile_options (-Wno-deprecated-builtins)
    add_compile_options (-Wno-deprecated-copy)
    if (ENABLE_OpenMP)
        if (CMAKE_HOST_APPLE)
            message(FATAL_ERROR "Apple Clang does not support OpenMP!")
        else ()
            add_compile_options (-fopenmp)
            add_link_options (-fopenmp)
        endif ()
    endif ()

elseif (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    # An error in Trilinos/12.18.1 causes
    # -Waggressive-loop-optimizations with -O3 optimization.
    # See issue 587 and corresponding merge request 415.
    # Clang does not support this flag!
    add_compile_options (-fno-aggressive-loop-optimizations)

    if (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.7.3")
        message (FATAL_ERROR "To build OPAL you need GCC version 4.7.3 or greater")
    endif ()
    # Warning flags

    if (ENABLE_OpenMP)
        add_compile_options (-fopenmp)
        add_link_options (-fopenmp)
    endif ()

    # Enables extra error checking in the form of precondition assertion, such
    # as bounds checking and null pointer checks when dereferencing smart pointers
    if (NOT ${BUILD_FOR_PRODUCTION})
        add_compile_options(-Wp,-D_GLIBCXX_ASSERTIONS)
    endif()
    enable_language (Fortran OPTIONAL)
else ()
    message (STATUS "Unknown C++ compiler. Please use the GNU or Intel compiler, if you are having problems.")
endif ()

message (STATUS "Fortran compilers: ${CMAKE_Fortran_COMPILER_NAMES}")
if (CMAKE_Fortran_COMPILER)
    message (STATUS "The Fortran compiler identification is: ${CMAKE_Fortran_COMPILER_ID}")
    message (STATUS "The Fortran compiler version is: ${CMAKE_Fortran_COMPILER_VERSION}")
    message (STATUS "The Fortran compiler is: ${CMAKE_Fortran_COMPILER}")
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        set (Fortran_LIBS "-lgfortran")
    endif ()
else ()
    message (STATUS "No Fortran compiler found (optional)")
    if (ENABLE_AMR)
        message (FATAL_ERROR "AMR is enabled but no Fortran compiler found!")
    endif ()
endif ()

find_package (H5Hut REQUIRED)

if (NOT ${WILL_BUILD_SHARED_LIBRARY})
    set (HDF5_USE_STATIC_LIBRARIES ON)
endif ()
set (HDF5_PREFER_PARALLEL ON)
find_package (HDF5 REQUIRED)
if(UNIX AND NOT APPLE)
    set (HDF5_LIBRARIES ${HDF5_LIBRARIES} ${CMAKE_DL_LIBS} rt z)
endif ()
message (STATUS "Found HDF5 libraries: ${HDF5_LIBRARIES}")

find_package (GSL REQUIRED)
message (STATUS "Found gsl libraries: ${GSL_LIBRARIES}")

set (BOOSTROOT $ENV{BOOST_DIR})
if (NOT ${WILL_BUILD_SHARED_LIBRARY})
    set (Boost_USE_STATIC_LIBS ON)
endif ()
set (Boost_USE_MULTITHREADED OFF)
set (Boost_USE_STATIC_RUNTIME OFF)
set (Boost_LIBS filesystem iostreams regex serialization system timer)
find_package (Boost 1.66.0
    REQUIRED COMPONENTS ${Boost_LIBS})
if(UNIX AND NOT APPLE)
    set (Boost_LIBRARIES ${Boost_LIBRARIES} rt)
endif ()
message (STATUS "Found Boost include dirs: ${Boost_INCLUDE_DIRS}")
message (STATUS "Found Boost library dirs: ${Boost_LIBRARY_DIRS}")
message (STATUS "Found Boost libraries: ${Boost_LIBRARIES}")
include_directories (SYSTEM ${Boost_INCLUDE_DIRS})

# Handle options
option (BUILD_OPAL_UNIT_TESTS "Unit tests" OFF)
if (BUILD_OPAL_UNIT_TESTS)
    find_package (GTest REQUIRED HINT $ENV{GTEST_PREFIX})
    add_definitions (-DWITH_UNIT_TESTS)
endif ()

option (ENABLE_AMR "Enable AMReX based AMR solver" OFF)
if (ENABLE_AMR)
    message (STATUS "Enable AMR_SOLVER    " ${ENABLE_AMR})
    enable_language (Fortran)

    find_package (AMReX CONFIG REQUIRED
        HINTS $ENV{AMREX_HOME}/cmake $ENV{AMREX_DIR}/cmake $ENV{AMREX_PREFIX}/cmake)

    message (STATUS "Found AMReX include dir: ${AMREX_INCLUDE_DIR}")
    message (STATUS "Found AMReX library dir: ${AMREX_LIBRARY_DIR}")
    message (STATUS "Found AMReX libraries: ${AMREX_LIBRARIES}")
    message (STATUS "Found AMReX Fortran links: ${AMREX_EXTRA_Fortran_LINK_LINE}")
    echo_amrex_config_options()

    add_definitions(${AMREX_DEFINES})
    add_definitions(-DENABLE_AMR)
endif ()

option (ENABLE_AMR_MG_SOLVER "Enable AMR MG solver" OFF)
option (ENABLE_SAAMG_SOLVER "Enable iterative SA-AMG-PCG self field solver" OFF)
if (ENABLE_SAAMG_SOLVER OR ENABLE_AMR_MG_SOLVER)
    #
    # :FIXME: We have to add checks whether the required Trilinos modules
    #         are available. If there are problems of any kind with this
    #         configuration, either fix it in a general way or build your
    #         Trilinos. Host specific should be avoided in almost any case.
    message (STATUS "Enable SAAMG_SOLVER: " ${ENABLE_SAAMG_SOLVER})
    message (STATUS "Enable AMR_MG_SOLVER: " ${ENABLE_AMR_MG_SOLVER})

    find_package (Trilinos REQUIRED HINTS $ENV{TRILINOS_PREFIX} $ENV{TRILINOS_DIR} $ENV{TRILINOS})

    message (STATUS "Found Trilinos: ${Trilinos_DIR}")
    message (STATUS "      Trilinos version: ${Trilinos_VERSION}")
    message (STATUS "      Trilinos package list: ${Trilinos_PACKAGE_LIST}")
    message (STATUS "      Trilinos libraries: ${Trilinos_LIBRARIES}")
    message (STATUS "      Trilinos TPL libraries: ${Trilinos_TPL_LIBRARIES}")

    if (NOT ${MPI_CXX_COMPILER} STREQUAL ${Trilinos_CXX_COMPILER} )
        message (WARNING "Compiler mismatch:")
        message (WARNING "        Trilinos was compiled with:")
        message (WARNING "            ${Trilinos_C_COMPILER}")
        message (WARNING "            ${Trilinos_CXX_COMPILER}")
        message (WARNING "            ${Trilinos_Fortran_COMPILER}")
        message (WARNING "        You are using:")
        message (WARNING "            ${CMAKE_C_COMPILER}")
        message (WARNING "            ${CMAKE_CXX_COMPILER}")
        message (WARNING "            ${CMAKE_Fortran_COMPILER}")
        message (FATAL_ERROR  "")
    endif ()
endif ()

option (ENABLE_OPAL_FEL "Enable OPAL FEL" OFF)
if (ENABLE_OPAL_FEL)
    message (STATUS "Enable OPAL FEL: " ${ENABLE_OPAL_FEL})
    find_package (MITHRA MODULE REQUIRED)
    add_definitions (-DENABLE_OPAL_FEL)
endif()

option (DBG_SCALARFIELD "Enable dump of scalar field rho_m" OFF)
if (DBG_SCALARFIELD)
    message (STATUS "Write scalar rho_m field is enabled ")
    add_definitions (-DDBG_SCALARFIELD)
endif ()

option (DBG_STENCIL "Dump stencil of iterative solver in matlab format" OFF)
if (DBG_STENCIL)
    message (STATUS "Write stencil is enabled ")
    add_definitions(-DDBG_STENCIL)
endif ()

option (NO_FIELD_ASSIGN_OPTIMIZATION "Disable compiler optimization of IPPL field assignment" OFF)
if (NO_FIELD_ASSIGN_OPTIMIZATION)
    message (STATUS "Compiler will not try to optimize field assignment")
    add_definitions(-DdontOPTIMIZE_FIELD_ASSIGNMENT)
endif ()

option (ENABLE_DOXYDOC "compile Doxygen documentation" OFF)
if (ENABLE_DOXYDOC)
    find_package(Doxygen REQUIRED)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY
	)
    add_custom_target (doxydoc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif ()

add_custom_target(always_rebuild ALL
    DEPENDS
    src/OPALrevision.h
    )

add_custom_command(
    OUTPUT src/OPALrevision.h
    COMMAND ./gen_OPALrevision
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS
    always_rebuild
    )

add_subdirectory (ippl)
add_subdirectory (optimizer)
add_subdirectory (src)
add_subdirectory (tools)

if (BUILD_OPAL_UNIT_TESTS)
    add_subdirectory (tests)
endif ()

# make variables visible in other CMakeLists files
set (OPAL_CXX_FLAGS ${OPAL_CXX_FLAGS}
    CACHE INTERNAL "" FORCE
    )
set (OPAL_LIBS ${OPAL_LIBS}
    CACHE INTERNAL "" FORCE
    )

# configure file
configure_file (
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake
    )

install (
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake
    DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/${PROJECT_NAME}"
    RENAME ${PROJECT_NAME}Config.cmake
    )