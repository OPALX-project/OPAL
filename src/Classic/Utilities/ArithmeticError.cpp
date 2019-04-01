// ------------------------------------------------------------------------
// $RCSfile: ArithmeticError.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ArithmeticError
//   The abstract base class for all CLASSIC arithmetic exceptions.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/ArithmeticError.h"


// Class ArithmeticError
// ------------------------------------------------------------------------


ArithmeticError::ArithmeticError
(const std::string &meth, const std::string &msg):
    ClassicException(meth, msg)
{}


ArithmeticError::ArithmeticError(const ArithmeticError &rhs):
    ClassicException(rhs)
{}


ArithmeticError::~ArithmeticError()
{}
