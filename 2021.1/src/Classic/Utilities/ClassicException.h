#ifndef CLASSIC_ClassicException_HH
#define CLASSIC_ClassicException_HH

// ------------------------------------------------------------------------
// $RCSfile: ClassicException.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ClassicException
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include <string>

// Class ClassicException
// ------------------------------------------------------------------------
/// The abstract base class for all exceptions in CLASSIC.

class ClassicException {

public:

    /// Return the message string for the exception.
    virtual const std::string &what() const;

    /// Return the name of the method or function which detected the exception.
    virtual const std::string &where() const;

protected:

    /// The usual constructor.
    // Arguments:
    // [DL]
    // [DT][b]meth[/b]
    // [DD]the name of the method or function detecting the exception
    // [DT][b]msg [/b]
    // [DD]the message string identifying the exception
    // [/DL]
    ClassicException(const std::string &meth, const std::string &msg);

    ClassicException(const ClassicException &);
    virtual ~ClassicException();

protected:

    // Not implemented.
    ClassicException();

    // The method detecting the exception and the message.
    const std::string message;
    const std::string method;
};

#endif // CLASSIC_ClassicException_HH
