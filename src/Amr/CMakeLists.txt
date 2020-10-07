set (_SRCS
    AmrBoxLib.cpp
    AmrObject.cpp
    AmrYtWriter.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

add_compile_options (-Wno-unused-variable
    -Wno-unused-but-set-variable
    -Wno-maybe-uninitialized
    )

set (HDRS
    AbstractAmrWriter.h
    AmrBoxLib.h
    AmrObject.h
    AmrYtWriter.h
    BoxLibLayout.h
    BoxLibLayout.hpp
    BoxLibParticle.h
    BoxLibParticle.hpp
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Amr")
