set (_SRCS
  Action.cpp
  Attribute.cpp
  AttributeBase.cpp
  AttributeHandler.cpp
  BeamSequence.cpp
  Definition.cpp
  Directory.cpp
  Editor.cpp
  Element.cpp
  Invalidator.cpp
  OpalData.cpp
  Object.cpp
  ObjectFunction.cpp
  PlaceRep.cpp
  RangeRep.cpp
  Table.cpp
  TableRowRep.cpp
  ValueDefinition.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_opal_sources (${_SRCS})

set (HDRS
    Action.h
    AttributeBase.h
    Attribute.h
    AttributeHandler.h
    BeamSequence.h
    Definition.h
    Directory.h
    Editor.h
    Element.h
    Expressions.h
    Invalidator.h
    ObjectFunction.h
    Object.h
    OpalData.h
    PlaceRep.h
    RangeRep.h
    Table.h
    TableRowRep.h
    ValueDefinition.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/AbstractObjects")

# vi: set et ts=4 sw=4 sts=4:
  
# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
