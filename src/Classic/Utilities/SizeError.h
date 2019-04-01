#ifndef CLASSIC_SizeError_HH
#define CLASSIC_SizeError_HH

// ------------------------------------------------------------------------
// $RCSfile: SizeError.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SizeError
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

// Class SizeError
// ------------------------------------------------------------------------
/// Size error exception.
//  This exception is thrown, when a CLASSIC routine or method is called
//  with arrays of inconsistent dimensions.

class SizeError: public ArithmeticError {

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
    SizeError(const std::string &meth, const std::string &msg);

    SizeError(const SizeError &);
    virtual ~SizeError();

private:

    // Not implemented.
    SizeError();
};

#endif // CLASSIC_SizeError_HH
