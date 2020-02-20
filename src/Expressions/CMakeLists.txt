set (_SRCS
    AColumn.cpp
    ARow.cpp
    ATable.cpp
    Expressions.cpp
    SCell.cpp
    SFunction.cpp
    SHash.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    AAutomatic.h
    ABinary.h
    AColumn.h
    ADeferred.h
    AList.h
    ARefAttr.h
    ARefExpr.h
    ARow.h
    ASUnary.h
    ATable.h
    AUnary.h
    AValue.h
    Expressions.h
    Indexer.h
    SAutomatic.h
    SBinary.h
    SCell.h
    SConstant.h
    SDeferred.h
    SFunction.h
    SHash.h
    SNull.h
    SRefAttr.h
    SRefExpr.h
    SUnary.h
    SValue.h
    TFind.h
    TFunction0.h
    TFunction1.h
    TFunction2.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Expressions")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
