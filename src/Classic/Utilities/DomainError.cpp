// ------------------------------------------------------------------------
// $RCSfile: DomainError.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DomainError
//   Argument(s) of mathematical function outside valid domain.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/DomainError.h"


// Class DomainError
// ------------------------------------------------------------------------


DomainError::DomainError(const std::string &meth):
    ArithmeticError(meth, "Domain error.")
{}


DomainError::DomainError(const DomainError &rhs):
    ArithmeticError(rhs)
{}


DomainError::~DomainError()
{}
