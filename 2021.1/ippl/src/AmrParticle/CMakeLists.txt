set (_SRCS
    )

set (_HDRS
    ParticleAmrLayout.h
    AmrParticleBase.h
    AmrParticleBase.hpp
    AmrParticleLevelCounter.h
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

ADD_IPPL_SOURCES (${_SRCS})
ADD_IPPL_HEADERS (${_HDRS})

install (FILES ${_HDRS} DESTINATION include/AmrParticle)
