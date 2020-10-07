set (_SRCS
  )

set (_HDRS
  BasicPairBuilder.h
  HashPairBuilder.h
  HashPairBuilderPeriodic.h
  PairConditions.h
  SortingPairBuilder.h
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_ippl_sources (${_SRCS})
add_ippl_headers (${_HDRS})

install (FILES ${_HDRS} DESTINATION include/Particle/PairBuilder)
