cmake_minimum_required (VERSION 3.1)
project (OPAL VERSION 2.1.0)
set (PROJECT_BUGREPORT opal@lists.psi.ch)

string (REGEX REPLACE "\\.([0-9])\\." ".0\\1." PACKAGE_VERSION ${PROJECT_VERSION})
string (REGEX REPLACE "\\.([0-9])$" ".0\\1" PACKAGE_VERSION ${PACKAGE_VERSION})
string (REGEX REPLACE "\\." "" PACKAGE_VERSION ${PACKAGE_VERSION})

message (PACKAGE_VERSION=${PACKAGE_VERSION})

if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
	)
endif ()

cmake_host_system_information (RESULT HOSTNAME QUERY HOSTNAME)

message ("Build type is: " ${CMAKE_BUILD_TYPE})
message ("Host OS System: ${CMAKE_HOST_SYSTEM}")
message ("Hostname:  ${HOSTNAME}")

# Select flags.
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")
set (CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

# Resolve all library dependencies
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

find_package (MPI REQUIRED)
message ("The C++ compiler identification is: ${CMAKE_CXX_COMPILER_ID}")
message ("The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
message ("The MPI C++ compiler is: ${MPI_CXX_COMPILER}")
message ("The underlying C++ compiler is: ${CMAKE_CXX_COMPILER}")

option (ENABLE_OpenMP "Use hybrid parallelism MPI-OpenMP" OFF)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

add_definitions (-DIPPL_LINUX -DIPPL_STRINGSTREAM)
add_definitions (-DIPPL_MPI -DMPICH_SKIP_MPICXX -DIPPL_DONT_POOL)
add_definitions (-DIPPL_USE_XDIV_RNG -DPETE_BITWISE_COPY)
add_definitions (-DIPPL_USE_PARTIAL_SPECIALIZATION)
add_definitions (-Drestrict=__restrict__ -DNOCTAssert)

add_compile_options (-Wall)
add_compile_options (-funroll-loops -fstrict-aliasing)
option (USE_STATIC_LIBRARIES "Link with static libraries if available" ON)
if (NOT ${USE_STATIC_LIBRARIES})
    add_compile_options (-fPIE -fPIC)
    set (CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -pie")
endif ()

# compiler dependent flags
message ("The underlying C++ compiler ID is: ${CMAKE_CXX_COMPILER_ID}")
if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel" )
    add_compile_options ("-diag-disable 383")
    add_compile_options ("-diag-disable 981")
    if (ENABLE_OpenMP)
	add_compile_options (-openmp)
    endif ()

elseif (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    add_compile_options (-Werror)
    add_compile_options (-Wsign-compare)
    add_compile_options (-Wunused-variable)
    add_compile_options (-Warray-bounds)
    # Increase template depth for boost
    add_compile_options (-ftemplate-depth=1024)
    # Allow overloaded virtual functions (instances to be fixed in OPAL)
    #add_compile_options (-Wno-overloaded-virtual)
    add_definitions(-DIPPL_RESTRICT_BUG)
    if (ENABLE_OpenMP)
	    add_compile_options (-fopenmp)
    endif ()

elseif (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    if (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.7.3")
        message (FATAL_ERROR "To build OPAL you need GCC version 4.7.3 or greater")
    endif ()
    # Warning flags, -Wno-unused-local-typedefs is needed for boost versions below 1.59
    add_compile_options (-Wno-unused-local-typedefs -Werror=unused-variable)

    if (ENABLE_OpenMP)
        add_compile_options (-fopenmp)
    endif ()

else ()
    message ("Unknown C++ compiler. Please use the GNU or Intel compiler, if you are having problems.")
endif ()

enable_language (Fortran OPTIONAL)
message ("Fortran compilers: ${CMAKE_Fortran_COMPILER_NAMES}")
if (CMAKE_Fortran_COMPILER)
    message ("The Fortran compiler identification is: ${CMAKE_Fortran_COMPILER_ID}")
    message ("The Fortran compiler version is: ${CMAKE_Fortran_COMPILER_VERSION}")
    message ("The Fortran compiler is: ${CMAKE_Fortran_COMPILER}")
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        set (Fortran_LIBS "-lgfortran")
    endif ()
else ()
    message ("No Fortran compiler found (optional)")
    if (ENABLE_AMR)
	message (FATAL_ERROR "AMR is enabled but no Fortran compiler found!")
    endif ()
endif ()

find_package (H5Hut REQUIRED)

if (USE_STATIC_LIBRARIES)
    set (HDF5_USE_STATIC_LIBRARIES ON)
endif ()
find_package (HDF5 REQUIRED)
message ("Found HDF5 libraries: ${HDF5_C_LIBRARIES}")

find_package (GSL REQUIRED)
message ("Found gsl libraries: ${GSL_LIBRARIES}")

set (BOOSTROOT $ENV{BOOST_DIR})
if (USE_STATIC_LIBRARIES)
    set (Boost_USE_STATIC_LIBS ON)
endif ()
set (Boost_USE_MULTITHREADED OFF)
set (Boost_USE_STATIC_RUNTIME OFF)
find_package (Boost 1.60.0
    REQUIRED COMPONENTS chrono filesystem iostreams regex serialization system timer)
message ("Found Boost include dir: ${Boost_INCLUDE_DIR}")
message ("Found Boost library dir: ${Boost_LIBRARY_DIR}")
message ("Found Boost libraries: ${Boost_LIBRARIES}")
include_directories (${Boost_INCLUDE_DIRS})

option (ENABLE_DKS "Enable DKS" OFF)
if (ENABLE_DKS)
    find_package (DKS 1.1.1 REQUIRED HINTS $ENV{DKS_PREFIX} $ENV{DKS_DIR} $ENV{DKS})
    set (DKS_VERSION_OPAL \"${DKS_VERSION}\")
    message ("Found DKS version: ${DKS_VERSION}")
    message ("Found DKS library: ${DKS_LIBRARY}")
    message ("Found DKS include dir: ${DKS_INCLUDE_DIR}")
    add_definitions(-DDKS_MPI)

    ### OpenCL compiler flags ###
    #add_compile_options (-lOpenCL -pthread -DDKS_OPENCL)

    ### CUDA compiler flags ###
    add_definitions(-DDKS_CUDA)

    ### if any accelerator enabled set flag to use DKS ###
    add_definitions(-DIPPL_DKS -DIPPL_DKS_CUDA -DOPAL_DKS)

endif ()

# Handle options
option (BUILD_OPAL_UNIT_TESTS "Unit tests" OFF)
if (BUILD_OPAL_UNIT_TESTS)
    find_package (GTest REQUIRED HINT $ENV{GTEST_PREFIX})
    add_definitions (-DWITH_UNIT_TESTS)
endif ()

option (ENABLE_AMR "Enable AMReX based AMR solver" OFF)
if (ENABLE_AMR)
    message ("Enable AMR_SOLVER    " ${ENABLE_AMR})
    enable_language (Fortran)

    find_package (AMReX CONFIG REQUIRED
        HINTS $ENV{AMREX_HOME}/cmake $ENV{AMREX_DIR}/cmake $ENV{AMREX_PREFIX}/cmake)

    message ("Found AMReX include dir: ${AMREX_INCLUDE_DIR}")
    message ("Found AMReX library dir: ${AMREX_LIBRARY_DIR}")
    message ("Found AMReX libraries: ${AMREX_LIBRARIES}")
    message ("Found AMReX Fortran links: ${AMREX_EXTRA_Fortran_LINK_LINE}")
    echo_amrex_config_options()

    add_definitions(${AMREX_DEFINES})
    add_definitions(-DENABLE_AMR)
endif ()

option (ENABLE_AMR_MG_SOLVER "Enable AMR MG solver" OFF)
option (ENABLE_SAAMG_SOLVER "Enable iteartive SA-AMG-PCG self field solver" OFF)
if (ENABLE_SAAMG_SOLVER OR ENABLE_AMR_MG_SOLVER)
    #
    # :FIXME: We have to add checks whether the required Trilinos modules
    #         are available. If there are problems of any kind with this
    #         configuration, either fix it in a general way or build your
    #         Trilinos. Host specific should be avoided in almost any case.
    message ("Enable SAAMG_SOLVER: " ${ENABLE_SAAMG_SOLVER})
    message ("Enable AMR_MG_SOLVER: " ${ENABLE_AMR_MG_SOLVER})

    find_package (Trilinos REQUIRED HINTS $ENV{TRILINOS_PREFIX} $ENV{TRILINOS_DIR} $ENV{TRILINOS})

    message ("Found Trilinos: ${Trilinos_DIR}")
    message ("      Trilinos version: ${Trilinos_VERSION}")
    message ("      Trilinos package list: ${Trilinos_PACKAGE_LIST}")
    message ("      Trilinos libraries: ${Trilinos_LIBRARIES}")
    message ("      Trilinos TPL libraries: ${Trilinos_TPL_LIBRARIES}")

    if (NOT ${MPI_CXX_COMPILER} STREQUAL ${Trilinos_CXX_COMPILER} )
        message ("Compiler mismatch:")
        message ("        Trilinos was compiled with:")
        message ("            ${Trilinos_C_COMPILER}")
        message ("            ${Trilinos_CXX_COMPILER}")
        message ("            ${Trilinos_Fortran_COMPILER}")
        message ("        You are using:")
        message ("            ${CMAKE_C_COMPILER}")
        message ("            ${CMAKE_CXX_COMPILER}")
        message ("            ${CMAKE_Fortran_COMPILER}")
        message (FATAL_ERROR  "")
    endif ()
endif ()

option (ENABLE_OPAL_FEL "Enable OPAL FEL" OFF)
if (ENABLE_OPAL_FEL)
    message ("Enable OPAL FEL: " ${ENABLE_OPAL_FEL})
    find_package (MITHRA MODULE REQUIRED)
    add_DEFINITIONS (-DOPAL_FEL)
endif()

option (DBG_SCALARFIELD "Enable dump of scalar field rho_m" OFF)
if (DBG_SCALARFIELD)
    message ("\nWrite scalar rho_m field is enabled ")
    add_DEFINITIONS (-DDBG_SCALARFIELD)
endif ()

option (DBG_STENCIL "Dump stencil of iterative solver in matlab format" OFF)
if (DBG_STENCIL)
    MESSAGE ("\nWrite stencil is enabled ")
    ADD_DEFINITIONS(-DDBG_STENCIL)
endif ()

option (NO_FIELD_ASSIGN_OPTIMIZATION "Disable compiler optimization of IPPL field assignment" OFF)
if (NO_FIELD_ASSIGN_OPTIMIZATION)
    MESSAGE ("\nCompiler will not try to optimize field assignment")
    ADD_DEFINITIONS(-DdontOPTIMIZE_FIELD_ASSIGNMENT)
endif ()

option (ENABLE_DOXYDOC "compile Doxygen documentation" OFF)

include (FindGit)
if (GIT_FOUND)
    set (OPAL_WC_REVISION 0)
    execute_PROCESS(
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMAND git rev-parse HEAD
        OUTPUT_VARIABLE GIT_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    set (OPAL_WC_REVISION \"${GIT_OUT}\")
else ()
    message ("could not find git")
    set (OPAL_WC_REVISION \"\")
endif ()

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

# vi: set et ts=4 sw=4 sts=4:
# Local Variables:
# mode:cmake
# cmake-tab-width: 4
# indent-tabs-mode:nil
# End:
