// ------------------------------------------------------------------------
// $RCSfile: AttributeError.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttributeError
//   The class for all CLASSIC exceptions related to object attributes.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Utilities/AttributeError.h"


// Class AttributeError
// ------------------------------------------------------------------------


AttributeError::AttributeError
(const std::string &meth, const std::string &msg):
    ClassicException(meth, msg)
{}


AttributeError::AttributeError(const AttributeError &rhs):
    ClassicException(rhs)
{}


AttributeError::~AttributeError()
{}
