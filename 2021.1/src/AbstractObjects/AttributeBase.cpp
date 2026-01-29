// ------------------------------------------------------------------------
// $RCSfile: AttributeBase.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttributeBase
//   The abstract base class for all attribute types.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:07 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/AttributeBase.h"
#include <sstream>

// Class AttributeBase
// ------------------------------------------------------------------------

AttributeBase::AttributeBase():
    RCObject()
{}


AttributeBase::~AttributeBase()
{}


std::string AttributeBase::getImage() const {
    std::ostringstream os;
    print(os);
    os << std::ends;
    return os.str();
}


bool AttributeBase::isExpression() const {
    return false;
}
