// ------------------------------------------------------------------------
// $RCSfile: String.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: String
//   A class used to parse string attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/String.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/SValue.h"
#include "Utilities/OpalException.h"
#include "Utilities/ParseError.h"

using namespace Expressions;


// Class String
// ------------------------------------------------------------------------

namespace Attributes {

    String::String(const std::string &name, const std::string &help):
        AttributeHandler(name, help, 0)
    {}


    String::~String()
    {}


    const std::string &String::getType() const {
        static const std::string type("string");
        return type;
    }


    void String::parse(Attribute &attr, Statement &stat, bool) const {

        std::string result = parseStringValue(stat, "String value expected.");
        attr.set(new SValue<std::string>(result));
    }

};