// ------------------------------------------------------------------------
// $RCSfile: ClassicException.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ClassicException
//   The base class for all CLASSIC exceptions.
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


// Class ClassicException
// ------------------------------------------------------------------------

ClassicException::ClassicException
(const std::string &meth, const std::string &msg):
    message(msg), method(meth)
{}


ClassicException::ClassicException(const ClassicException &rhs):
    message(rhs.message), method(rhs.method)
{}


ClassicException::~ClassicException()
{}


const std::string &ClassicException::what() const {
    return message;
}


const std::string &ClassicException::where() const {
    return method;
}
