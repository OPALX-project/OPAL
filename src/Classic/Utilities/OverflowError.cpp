// ------------------------------------------------------------------------
// $RCSfile: OverflowError.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OverflowError
//   Overflow in mathematical function.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/OverflowError.h"


// Class OverflowError
// ------------------------------------------------------------------------


OverflowError::OverflowError(const std::string &meth):
    ArithmeticError(meth, "Overflow error.")
{}


OverflowError::OverflowError(const OverflowError &rhs):
    ArithmeticError(rhs)
{}


OverflowError::~OverflowError()
{}
