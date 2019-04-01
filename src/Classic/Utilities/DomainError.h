#ifndef CLASSIC_DomainError_HH
#define CLASSIC_DomainError_HH

// ------------------------------------------------------------------------
// $RCSfile: DomainError.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DomainError
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

// Class DomainError
// ------------------------------------------------------------------------
/// Domain error exception.
//  This exception is thrown, when a domain error is detected in a function.

class DomainError: public ArithmeticError {

public:

    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [/DL]
    DomainError(const std::string &meth);

    DomainError(const DomainError &);
    virtual ~DomainError();

private:

    // Not implemented.
    DomainError();
};

#endif // CLASSIC_DomainError_HH
