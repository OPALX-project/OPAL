#
# Find MITHRA package
# https://github.com/aryafallahi/mithra
#
# MITHRA_INCLUDE_DIR
# MITHRA_LIBRARY_DIR
# MITHRA_FOUND

find_path (MITHRA_INCLUDE_DIR mithra/classes.h
    HINTS $ENV{MITHRA_INCLUDE_DIR} $ENV{MITHRA_INCLUDE_PATH} $ENV{MITHRA_PREFIX}/include $ENV{MITHRA}/include
)

find_path (MITHRA_LIBRARY_DIR libmithra.a
    HINTS $ENV{MITHRA_LIBRARY_DIR} $ENV{MITHRA_LIBRARY_PATH} $ENV{MITHRA_PREFIX}/lib $ENV{MITHRA}/lib
)

if (MITHRA_INCLUDE_DIR AND MITHRA_LIBRARY_DIR)
    set (MITHRA_FOUND "YES")
endif ()

if (MITHRA_FOUND)
    if (NOT MITHRA_FIND_QUIETLY)
        message (STATUS "Found MITHRA include dir: ${MITHRA_INCLUDE_DIR}")
        message (STATUS "Found MITHRA library dir: ${MITHRA_LIBRARY_DIR}")
    endif ()
else (MITHRA_FOUND)
    if (MITHRA_FIND_REQUIRED)
        if (NOT MITHRA_INCLUDE_DIR)
            message (WARNING 
                "MITHRA include directory was not found! "
                "Make sure that MITHRA is compiled and that "
                "the directory mithra/include/mithra has been automatically created. "
                "Also make sure that at least one of the following "
                "environment variables is set: "
                "MITHRA_INCLUDE_DIR, MITHRA_INCLUDE_PATH, MITHRA_PREFIX, or MITHRA.")
        endif ()
        if (NOT MITHRA_LIBRARY_DIR)
            message (WARNING 
                "MITHRA library was not found! "
                "Make sure that MITHRA is compiled and that "
                "the directory mithra/lib has been automatically created. "
                "Also make sure that at least one of the following "
                "environment variables is set: "
                "MITHRA_LIBRARY_DIR, MITHRA_LIBRARY_PATH, MITHRA_PREFIX, or MITHRA.")
        endif ()
        message (STATUS "MITHRA can be downloaded and compiled from https://github.com/aryafallahi/mithra.git")
        message (FATAL_ERROR "Could not find MITHRA!")
    endif (MITHRA_FIND_REQUIRED)
endif (MITHRA_FOUND)
