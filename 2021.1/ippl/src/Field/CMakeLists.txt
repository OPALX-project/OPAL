set (_SRCS
    )

set (_HDRS
    Assign.hpp
    Assign.h
    AssignDefs.h
    AssignGeneralBF.hpp
    AssignGeneralIBF.hpp
    AssignTags.h
    BareField.hpp
    BareField.h
    BareFieldIterator.h
    BCond.hpp
    BCond.h
    BrickExpression.hpp
    BrickExpression.h
    BrickIterator.hpp
    BrickIterator.h
    CompressedBrickIterator.hpp
    CompressedBrickIterator.h
    Field.hpp
    Field.h
    FieldLoc.h
    FieldSpec.h
    GuardCellSizes.hpp
    GuardCellSizes.h
    IndexedBareField.hpp
    IndexedBareField.h
    IndexedField.h
    LField.hpp
    LField.h
    ReductionLoc.hpp
    ReductionLoc.h
    )

include_DIRECTORIES (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_ippl_sources (${_SRCS})
add_ippl_headers (${_HDRS})

install (FILES ${_HDRS} DESTINATION include/Field)
