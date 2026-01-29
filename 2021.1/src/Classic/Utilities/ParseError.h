#ifndef CLASSIC_ParseError_HH
#define CLASSIC_ParseError_HH

// ------------------------------------------------------------------------
// $RCSfile: ParseError.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParseError
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


// Class ParseError
// ------------------------------------------------------------------------
/// Parse exception.
//  This exception is thrown by the CLASSIC parser whein it detects an
//  input format error.

class ParseError: public ClassicException {

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
    ParseError(const std::string &meth, const std::string &msg);

    ParseError(const ParseError &);
    virtual ~ParseError();

private:

    // Not implemented.
    ParseError();
};

#endif // CLASSIC_ParseError_HH
