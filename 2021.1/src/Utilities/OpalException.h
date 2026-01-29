#ifndef OPAL_OpalException_HH
#define OPAL_OpalException_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalException.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalException
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:48 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Utilities/ClassicException.h"


// Class OpalException
// ------------------------------------------------------------------------
/// The base class for all OPAL exceptions.

class OpalException: public ClassicException {

public:

    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [DT][b]msg [/b]
    // [DD]the message string identifying the exception
    // [/DL]
    explicit OpalException(const std::string &meth, const std::string &msg);

    OpalException(const OpalException &);
    virtual ~OpalException();

    /// Return the message string for the exception.
    using ClassicException::what;

    /// Return the name of the method or function which detected the exception.
    using ClassicException::where;

private:

    // Not implemented.
    OpalException();
};

#endif // OPAL_OpalException_HH
