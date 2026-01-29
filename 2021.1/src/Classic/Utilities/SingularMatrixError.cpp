// ------------------------------------------------------------------------
// $RCSfile: SingularMatrixError.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SingularMatrixError
//   Singular matrix in matrix inversion.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/SingularMatrixError.h"


// Class SingularMatrixError
// ------------------------------------------------------------------------


SingularMatrixError::SingularMatrixError(const std::string &meth):
    ArithmeticError(meth, "Singular matrix.")
{}


SingularMatrixError::SingularMatrixError(const SingularMatrixError &rhs):
    ArithmeticError(rhs)
{}


SingularMatrixError::~SingularMatrixError()
{}
