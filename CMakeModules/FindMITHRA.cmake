#
# Find MITHRA package
# https://github.com/aryafallahi/mithra
#
# MITHRA_INCLUDE_DIR
# MITHRA_FOUND

find_path (MITHRA_INCLUDE_DIR src/classes.h)

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
