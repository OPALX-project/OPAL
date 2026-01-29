// ------------------------------------------------------------------------
// $RCSfile: EigenvalueError.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: EigenvalueError
//   Unable to find all eigenvalues and/or eigenvectors.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/EigenvalueError.h"


// Class EigenvalueError
// ------------------------------------------------------------------------


EigenvalueError::EigenvalueError(const std::string &meth, const std::string &msg):
    ArithmeticError(meth, msg)
{}


EigenvalueError::EigenvalueError(const EigenvalueError &rhs):
    ArithmeticError(rhs)
{}


EigenvalueError::~EigenvalueError()
{}
