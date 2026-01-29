#
# Find GSL includes and libraries
#
# The FindGSL module shipped with CMake has some drawbacks. To be 
# able to find GSL, the CMake variable GSL_ROOT_DIR *must* be set.
# The environment variable C_INCLUDE_PATH and LIBRARY_PATH are
# *not* taken into account.
#
# The following variables will be set if GSL is found:
#
# GSL_INCLUDE_DIR	where to find GSL include files
# GSL_LIBRARY		GSL library to link against.
# GSL_CBLAS_LIBRARY 	GSL CBlas library to link against
# GSL_LIBRARIES		GSL libraries required for linking
# GSL_FOUND		set to True if GSL was found
#

if( DEFINED ENV{GSL_ROOT_DIR} )
  set( GSL_ROOT_DIR $ENV{GSL_ROOT_DIR} )
elseif( DEFINED ENV{GSL_DIR} )
  set( GSL_ROOT_DIR $ENV{GSL_DIR} )
elseif( DEFINED ENV{GSL_HOME} )
  set( GSL_ROOT_DIR $ENV{GSL_HOME} )
elseif( DEFINED ENV{GSL_PREFIX} )
  set( GSL_ROOT_DIR $ENV{GSL_PREFIX} )
elseif( DEFINED ENV{GSL} )
  set( GSL_ROOT_DIR $ENV{GSL} )
else()
  set( GSL_ROOT_DIR "/usr" )
endif()

find_path( GSL_INCLUDE_DIR gsl/gsl_fft.h
  HINTS ${GSL_ROOT_DIR}/include $ENV{GSL_INCLUDE_PATH} $ENV{GSL_INCLUDE_DIR}
  PATHS ENV C_INCLUDE_PATH
  )

find_library( GSL_LIBRARY gsl
  HINTS ${GSL_ROOT_DIR}/lib $ENV{GSL_LIBRARY_PATH} $ENV{GSL_LIBRARY_DIR}
  PATHS ENV LIBRARY_PATH
  )

find_library (GSL_CBLAS_LIBRARY gslcblas
  HINTS ${GSL_ROOT_DIR}/lib $ENV{GSL_LIBRARY_PATH} $ENV{GSL_LIBRARY_DIR}
  PATHS ENV LIBRARY_PATH
  )

if( GSL_INCLUDE_DIR AND GSL_LIBRARY AND GSL_CBLAS_LIBRARY)
  set( GSL_FOUND "YES" )
  set( GSL_LIBRARIES ${GSL_LIBRARY},${GSL_CBLAS_LIBRARY} )
endif()

if( GSL_FOUND )
  if( NOT GSL_FIND_QUIETLY )
    message( STATUS "Found GSL libraries: ${GSL_LIBRARY}")
    message( STATUS "Found GSL include dir: ${GSL_INCLUDE_DIR}")
  endif()
else()
  if( GSL_FIND_REQUIRED )
    message( FATAL_ERROR "Could not find GSL!" )
  endif()
endif()
