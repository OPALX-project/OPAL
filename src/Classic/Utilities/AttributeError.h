#ifndef CLASSIC_AttributeError_HH
#define CLASSIC_AttributeError_HH

// ------------------------------------------------------------------------
// $RCSfile: AttributeError.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttributeError
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/ClassicException.h"

#include <string>

// Class AttributeError
// ------------------------------------------------------------------------
/// Exception for bad object attribute.
//  This exception is thrown when bad object attributes are detected.

class AttributeError: public ClassicException {

public:

    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [DT][b]msg [/b]
    // [DD]the message string identifying the exception
    // [/DL]
    AttributeError(const std::string &meth, const std::string &msg);

    AttributeError(const AttributeError &);
    virtual ~AttributeError();

private:

    // Not implemented.
    AttributeError();
};

#endif // CLASSIC_AttributeError_HH
