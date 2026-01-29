// ------------------------------------------------------------------------
// $RCSfile: DivideError.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DivideError
//   Divide by zero error.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/DivideError.h"


// Class DivideError
// ------------------------------------------------------------------------


DivideError::DivideError(const std::string &method):
    ArithmeticError(method, "Divide error.")
{}


DivideError::DivideError(const DivideError &rhs):
    ArithmeticError(rhs)
{}


DivideError::~DivideError()
{}
