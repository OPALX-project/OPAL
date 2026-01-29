#ifndef CLASSIC_SingularMatrixError_HH
#define CLASSIC_SingularMatrixError_HH

// ------------------------------------------------------------------------
// $RCSfile: SingularMatrixError.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SingularMatrixError
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/ArithmeticError.h"

#include <string>

// Class SingularMatrixError
// ------------------------------------------------------------------------
/// Singular matrix exception.
//  This exception is thrown, when an attempt is made to invert a singular
//  matrix.

class SingularMatrixError: public ArithmeticError {

public:

    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [/DL]
    // Construction/destruction.
    explicit SingularMatrixError(const std::string &meth);

    SingularMatrixError(const SingularMatrixError &);
    virtual ~SingularMatrixError();

private:

    // Not implemented.
    SingularMatrixError();
};

#endif // CLASSIC_SingularMatrixError_HH
