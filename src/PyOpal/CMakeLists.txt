list(APPEND Boost_LIBS python37)
find_package (Boost
    REQUIRED COMPONENTS ${Boost_LIBS})

set (PYUTILS_SRCS Globals.cpp ExceptionTranslation.cpp)

# FindPython is supposed to set Python_SITEARCH and Python_SITELIB - but as far
# as I can tell they are empty; so try setting PYTHON_SITE_PACKAGES instead
# this will install to the system default python site-packages dir
# ignores CMAKE_INSTALL_DIR (but OPAL seems to not comply with that anyway).
# make install seems to be pretty broken for OPAL in general
#
# To get this to work, I had to do make install DESTDIR=path/to/doom
# and then replace certain directories with soft link to the actual directory
# where I wanted to install. Darkness and evil.
#

execute_process ( 
COMMAND
python -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE
)



function(python_module MODULE_NAME PYTHON_SRCS)
    MESSAGE(STATUS "Building python module ${MODULE_NAME} with sources ${PYTHON_SRCS}")
    MESSAGE(STATUS "and install target ${PYTHON_SITE_PACKAGES}/pyopal/")
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
    # could do something using setuptools instead? e.g.
    # https://stackoverflow.com/questions/42585210/extending-setuptools-extension-to-use-cmake-in-setup-py   
    install (TARGETS ${MODULE_NAME} DESTINATION "${PYTHON_SITE_PACKAGES}/pyopal/")
endfunction(python_module)

# for Parser, we call OpalMain; so must not link to Globals.cpp as we get
# multiple definitions of Opal global objects (i.e. ippl and gmsg)
set (PARSER_SRCS PyParser.cpp ExceptionTranslation.cpp)
python_module(parser "${PARSER_SRCS}")

set (FIELD_SRCS PyField.cpp ${PYUTILS_SRCS})
python_module(field "${FIELD_SRCS}")
