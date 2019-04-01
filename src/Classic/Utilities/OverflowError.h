#ifndef CLASSIC_OverflowError_HH
#define CLASSIC_OverflowError_HH

// ------------------------------------------------------------------------
// $RCSfile: OverflowError.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OverflowError
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

// Class OverflowError
// ------------------------------------------------------------------------
/// Overflow exception.
//  This exception is thrown, when a function detects an arithmetic overflow.

class OverflowError: public ArithmeticError {

public:

    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [/DL]
    // Construction/destruction.
    OverflowError(const std::string &meth);

    OverflowError(const OverflowError &);
    virtual ~OverflowError();

private:

    // Not implemented.
    OverflowError();
};

#endif // CLASSIC_OverflowError_HH
