set (_SRCS
  )

set (_HDRS
  SubParticleAssign.hpp
  SubParticleAssignDefs.h
  SubParticleAssign.h
  SubParticleAttrib.hpp
  SubParticleAttrib.h
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_ippl_sources (${_SRCS})
add_ippl_headers (${_HDRS})

install (FILES ${_HDRS} DESTINATION include/SubParticle)
