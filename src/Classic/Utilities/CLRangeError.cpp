// ------------------------------------------------------------------------
// $RCSfile: CLRangeError.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CLRangeError
//   Index out of range in matrix or vector operation.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/CLRangeError.h"


// Class CLRangeError
// ------------------------------------------------------------------------


CLRangeError::CLRangeError(const std::string &meth, const std::string &msg):
    ArithmeticError(meth, msg)
{}


CLRangeError::CLRangeError(const CLRangeError &rhs):
    ArithmeticError(rhs)
{}


CLRangeError::~CLRangeError()
{}
