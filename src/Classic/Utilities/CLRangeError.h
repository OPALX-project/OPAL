#ifndef CLASSIC_CLRangeError_HH
#define CLASSIC_CLRangeError_HH

// ------------------------------------------------------------------------
// $RCSfile: CLRangeError.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CLRangeError
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

// Class CLRangeError
// ------------------------------------------------------------------------
/// Range error.
//  This exception is thrown, when a CLASSIC routine or method detects an
//  index out of range.

class CLRangeError: public ArithmeticError {

public:

    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [DT][b]msg [/b]
    // [DD]the message string identifying the exception
    // [/DL]
    // Construction/destruction.
    CLRangeError(const std::string &meth, const std::string &msg);

    CLRangeError(const CLRangeError &);
    virtual ~CLRangeError();

private:

    // Not implemented.
    CLRangeError();
};

#endif // CLASSIC_CLRangeError_HH
