// ------------------------------------------------------------------------
// $RCSfile: SizeError.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SizeError
//   Size by zero error.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/SizeError.h"


// Class SizeError
// ------------------------------------------------------------------------


SizeError::SizeError(const std::string &meth, const std::string &msg):
    ArithmeticError(meth, msg)
{}


SizeError::SizeError(const SizeError &rhs):
    ArithmeticError(rhs)
{}


SizeError::~SizeError()
{}
