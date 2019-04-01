// ------------------------------------------------------------------------
// $RCSfile: NormalFormError.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: NormalFormError
//   Index out of range in Algebra access.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/NormalFormError.h"


// Class NormalFormError
// ------------------------------------------------------------------------


NormalFormError::NormalFormError(const std::string &meth,
                                 const std::string &msg):
    ArithmeticError(meth, msg)
{}


NormalFormError::NormalFormError(const NormalFormError &rhs):
    ArithmeticError(rhs)
{}


NormalFormError::~NormalFormError()
{}
