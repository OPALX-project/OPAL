list(APPEND Boost_LIBS python37)
find_package (Boost
    REQUIRED COMPONENTS ${Boost_LIBS})

set (PYUTILS_SRCS Globals.cpp ExceptionTranslation.cpp)


function(python_module MODULE_NAME PYTHON_SRCS)
    MESSAGE(STATUS "Building python module ${MODULE_NAME} with sources ${PYTHON_SRCS}")
    add_library(${MODULE_NAME} SHARED ${PYTHON_SRCS})
    set_target_properties (${MODULE_NAME} PROPERTIES PREFIX "")
    target_link_libraries( ${MODULE_NAME}
        OPAL
        ${OPTP_LIBS}
        ${OPTP_LIBRARY}
        ${IPPL_LIBRARY}
        ${GSL_LIBRARY}
        ${GSL_CBLAS_LIBRARY}
        ${H5Hut_LIBRARY}
        ${HDF5_LIBRARIES}
        ${Boost_LIBRARIES}
        ${Python3_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        m
        z
        "-rdynamic"
    )
    target_include_directories(${MODULE_NAME} PUBLIC
        ${Python3_INCLUDE_DIRS}
        ${IPPL_INCLUDE_DIR}
        ${H5Hut_INCLUDE_DIR}
        ${HDF5_INCLUDE_DIR}
        ${GSL_INCLUDE_DIR}
        ${OPAL_SOURCE_DIR}/src
        ${OPAL_SOURCE_DIR}/src/Classic
        ${OPAL_SOURCE_DIR}/optimizer
    )
    target_link_directories(${MODULE_NAME} PUBLIC ${OPAL_SOURCE_DIR}/src)
    install (TARGETS ${MODULE_NAME} DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
endfunction(python_module)

# for Parser, we call OpalMain; so must not link to Globals.cpp as we get
# multiple definitions of Opal global objects (i.e. ippl and gmsg)
set (PARSER_SRCS PyParser.cpp ExceptionTranslation.cpp)
python_module(parser "${PARSER_SRCS}")

set (FIELD_SRCS PyField.cpp ${PYUTILS_SRCS})
python_module(field "${FIELD_SRCS}")
