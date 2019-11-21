set (_SRCS
    ArithmeticError.cpp
    AttributeError.cpp
    CLRangeError.cpp
    ClassicException.cpp
    ClassicField.cpp
    ClassicRandom.cpp
    ComplexErrorFun.cpp
    ConvergenceError.cpp
    DivideError.cpp
    DomainError.cpp
    EigenvalueError.cpp
    FormatError.cpp
    GeneralClassicException.cpp
    InverseGauss.cpp
    LogicalError.cpp
    Mesher.cpp
    MSLang.cpp
    MSLang/ArgumentExtractor.cpp
    MSLang/BoundingBox.cpp
    MSLang/Difference.cpp
    MSLang/Ellipse.cpp
    MSLang/Intersection.cpp
    MSLang/Mask.cpp
    MSLang/matheval.cpp
    MSLang/Polygon.cpp
    MSLang/QuadTree.cpp
    MSLang/Rectangle.cpp
    MSLang/Repeat.cpp
    MSLang/Rotation.cpp
    MSLang/Shear.cpp
    MSLang/SymmetricDifference.cpp
    MSLang/Translation.cpp
    MSLang/Triangle.cpp
    MSLang/Union.cpp
    NormalFormError.cpp
    Options.cpp
    OverflowError.cpp
    ParseError.cpp
    PortableBitmapReader.cpp
    PortableGraymapReader.cpp
    RingSection.cpp
    SingularMatrixError.cpp
    SizeError.cpp
    SwitcherError.cpp
    Util.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    ArithmeticError.h
    AttributeError.h
    CLRangeError.h
    ClassicException.h
    ClassicField.h
    ClassicRandom.h
    ComplexErrorFun.h
    ConvergenceError.h
    DivideError.h
    DomainError.h
    EigenvalueError.h
    FormatError.h
    GeneralClassicException.h
    InverseGauss.h
    LogicalError.h
    Mesher.h
    MSLang.h
    MSLang/ArgumentExtractor.h
    MSLang/AffineTransformation.h
    MSLang/BoundingBox.h
    MSLang/Difference.h
    MSLang/Ellipse.h
    MSLang/Intersection.h
    MSLang/Mask.h
    MSLang/matheval.h
    MSLang/matheval.hpp
    MSLang/Polygon.h
    MSLang/QuadTree.h
    MSLang/Rectangle.h
    MSLang/Repeat.h
    MSLang/Rotation.h
    MSLang/Shear.h
    MSLang/SymmetricDifference.h
    MSLang/Translation.h
    MSLang/Triangle.h
    MSLang/Union.h
    NormalFormError.h
    OptionTypes.h
    Options.h
    OverflowError.h
    ParseError.h
    PortableBitmapReader.h
    PortableGraymapReader.h
    RingSection.h
    SingularMatrixError.h
    SizeError.h
    SwitcherError.h
    Util.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Utilities")

# vi: set et ts=4 sw=4 sts=4:

# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# End:
