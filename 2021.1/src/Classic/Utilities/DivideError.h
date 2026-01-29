#ifndef CLASSIC_DivideError_HH
#define CLASSIC_DivideError_HH

// ------------------------------------------------------------------------
// $RCSfile: DivideError.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DivideError
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/ArithmeticError.h"

#include <string>

// Class DivideError
// ------------------------------------------------------------------------
/// Zero divide error.
//  This exception is thrown, when a division by zero occurs.

class DivideError: public ArithmeticError {

public:

    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [/DL]
    explicit DivideError(const std::string &meth);

    DivideError(const DivideError &);
    virtual ~DivideError();

private:

    // Not implemented.
    DivideError();
};

#endif // CLASSIC_DivideError_HH
