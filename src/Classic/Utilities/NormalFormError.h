#ifndef CLASSIC_NormalFormError_HH
#define CLASSIC_NormalFormError_HH

// ------------------------------------------------------------------------
// $RCSfile: NormalFormError.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: NormalFormError
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

// Class NormalFormError
// ------------------------------------------------------------------------
/// Normal form error exception.
//  This exception is thrown, when a normal form routine fails to find the
//  complete normal form, or when it detects an inconsistent call.

class NormalFormError: public ArithmeticError {

public:

    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [DT][b]msg [/b]
    // [DD]the message string identifying the exception
    // [/DL]
    NormalFormError(const std::string &meth, const std::string &msg);

    NormalFormError(const NormalFormError &);
    virtual ~NormalFormError();

private:

    // Not implemented.
    NormalFormError();
};

#endif // CLASSIC_NormalFormError_HH
