set (_SRCS
    fftpack.cpp
    )

set (_HDRS
    FFTBase.hpp
    FFTBase.h
    FFT.hpp
    FFT.h
    fftpack.h
    fftpack_FFT.h
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_ippl_sources (${_SRCS})
add_ippl_headers (${_HDRS})
