#
# Find MITHRA package
# https://github.com/aryafallahi/mithra
#
# MITHRA_INCLUDE_DIR
# MITHRA_FOUND

find_path (MITHRA_INCLUDE_DIR mithra/classes.hh
    HINTS $ENV{MITHRA_INCLUDE_PATH} $ENV{MITHRA_INCLUDE_DIR} $ENV{MITHRA_PREFIX}/include $ENV{MITHRA_DIR}/include $ENV{MITHRA}/include
    PATHS ENV C_INCLUDE_PATH CPLUS_INCLUDE_PATH
)

if (MITHRA_INCLUDE_DIR)
    set (MITHRA_FOUND "YES")
endif ()

if (MITHRA_FOUND)
    if (NOT MITHRA_FIND_QUIETLY)
        message (STATUS "Found MITHRA include dir: ${MITHRA_INCLUDE_DIR}")
    endif ()
else (MITHRA_FOUND)
    if (MITHRA_FIND_REQUIRED)
        message (FATAL_ERROR "Could not find MITHRA!")
    endif (MITHRA_FIND_REQUIRED)
endif (MITHRA_FOUND)
