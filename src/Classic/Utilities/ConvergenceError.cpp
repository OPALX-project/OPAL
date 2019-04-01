// ------------------------------------------------------------------------
// $RCSfile: ConvergenceError.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConvergenceError
//   Problems with convergence of iterations and series.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/ConvergenceError.h"


// Class ConvergenceError
// ------------------------------------------------------------------------


ConvergenceError::ConvergenceError(const std::string &meth, const std::string &msg):
    ArithmeticError(meth, msg)
{}


ConvergenceError::ConvergenceError(const ConvergenceError &rhs):
    ArithmeticError(rhs)
{}


ConvergenceError::~ConvergenceError()
{}
