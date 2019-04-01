#ifndef CLASSIC_FormatError_HH
#define CLASSIC_FormatError_HH

// ------------------------------------------------------------------------
// $RCSfile: FormatError.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FormatError
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/ClassicException.h"

#include <string>

// Class FormatError
// ------------------------------------------------------------------------
/// Format error exception.
//  This exception is thrown, when an input routine detects a format error.

class FormatError: public ClassicException {

public:

    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [DT][b]msg [/b]
    // [DD]the message string identifying the exception
    // [/DL]
    FormatError(const std::string &meth, const std::string &msg);

    FormatError(const FormatError &);
    virtual ~FormatError();

private:

    // Not implemented.
    FormatError();
};

#endif // CLASSIC_FormatError_HH
