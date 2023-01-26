# We need to redo the dependency for boost to include python. Because we don't
# want opal main functions to depend on Python.h, we redo the boost find_package 
# here but including boost::python now

# boost::python has the python version hardcoded in the library name (yuck) so
# we resolve the version.

SET(PyBoostLib "python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}")
MESSAGE(STATUS "Adding ${PyBoostLib} to boost required list")

# pyopal will be installed in PYOPAL_INSTALL_DIR; can be defined at cmake time
# or we look for the output of cmake finder; or we search using a simple python
# script
if (PYOPAL_INSTALL_DIR)
    # do nothing
else() # elif doesn't exist in cmake 3.16
    if (Python3_SITEARCH)
        # In my cmake build Python3_SITEARCH was not set
        SET(PYOPAL_INSTALL_DIR ${Python3_SITEARCH})
    else()
        exec_program("python find_install_location.py" "${CMAKE_CURRENT_SOURCE_DIR}"
                     OUTPUT_VARIABLE Python_INSTALL_LOCATION)
        SET(PYOPAL_INSTALL_DIR ${Python_INSTALL_LOCATION})
    endif()
endif()
MESSAGE(STATUS "PyOpal install location resolved to ${PYOPAL_INSTALL_DIR}")

set (Boost_LIBS ${Boost_LIBS} ${PyBoostLib})
find_package (Boost 1.66.0
    REQUIRED COMPONENTS ${Boost_LIBS})

SET(PYTHON_ALLOWED_INSTALL_DIRS objects elements)

# This is a generic "Make a python module" function that can be called to make
# a python module from C++ code like 
#            python_module(<module_name> <list of source files> <install target>)
# where <install target> must be one of ${PYTHON_ALLOWED_INSTALL_DIRS}
function(python_module MODULE_NAME PYTHON_SRCS INSTALL_DIRECTORY)
    # PYTHON_ALLOWED_INSTALL_DIRS is the list of the python modules that exists
    # on python side i.e.
    #     pyopal.<blah>
    SET(IS_ALLOWED "CHEESE")
    LIST(FIND PYTHON_ALLOWED_INSTALL_DIRS ${INSTALL_DIRECTORY} IS_ALLOWED)
    if("${IS_ALLOWED}" STREQUAL "-1")
        MESSAGE(FATAL_ERROR "Python module install dir '${INSTALL_DIRECTORY}' must be one of '${PYTHON_ALLOWED_INSTALL_DIRS}'")
    endif()
    add_library(${MODULE_NAME} SHARED ${PYTHON_SRCS})
    set_target_properties (${MODULE_NAME} PROPERTIES PREFIX "")
    target_link_libraries( ${MODULE_NAME}
        libOPALdynamic
        ${Python3_LIBRARIES}
        ${OPTP_LIBS}
        ${OPTP_LIBRARY}
        ${IPPL_LIBRARY}
        ${GSL_LIBRARY}
        ${GSL_CBLAS_LIBRARY}
        ${H5Hut_LIBRARY}
        ${HDF5_LIBRARIES}
        ${Boost_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        m
        z
        "-rdynamic"
    )
    target_include_directories(${MODULE_NAME} PUBLIC ${Python3_INCLUDE_DIRS})
    SET(MY_INSTALL_LOCATION ${PYOPAL_INSTALL_DIR}/pyopal/${INSTALL_DIRECTORY})
    MESSAGE(STATUS "Building python module ${MODULE_NAME} with sources ${PYTHON_SRCS}")
    install (TARGETS ${MODULE_NAME} DESTINATION "${MY_INSTALL_LOCATION}")
endfunction(python_module)

# This is a generic "Make a python module" function that can be called to make
# a python module from python code like
#            python_module(<python file> <install target>)
# where <install target> must be one of ${PYTHON_ALLOWED_INSTALL_DIRS}.
function(python_py_module PYTHON_LIB INSTALL_DIRECTORY)
    SET(IS_ALLOWED "CHEESE")
    LIST(FIND PYTHON_ALLOWED_INSTALL_DIRS ${INSTALL_DIRECTORY} IS_ALLOWED)
    if("${IS_ALLOWED}" STREQUAL "-1")
        MESSAGE(FATAL_ERROR "Python module install dir '${INSTALL_DIRECTORY}' must be one of '${PYTHON_ALLOWED_INSTALL_DIRS}'")
    endif()
    SET(MY_INSTALL_LOCATION ${PYOPAL_INSTALL_DIR}/pyopal/${INSTALL_DIRECTORY})
    MESSAGE(STATUS "Building python module ${MODULE_NAME} from ${PYTHON_LIB}")
    install (PROGRAMS ${PYTHON_LIB} DESTINATION "${MY_INSTALL_LOCATION}")
endfunction(python_py_module)

set(OPAL_SRCS "${OPAL_SRCS}" PARENT_SCOPE)

add_subdirectory(PyCore)
add_subdirectory(PyElements)
add_subdirectory(PyObjects)
add_subdirectory(PyPython)
