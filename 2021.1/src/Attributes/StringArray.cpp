// ------------------------------------------------------------------------
// $RCSfile: StringArray.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class StringArray:
//   A class used to parse array attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/StringArray.h"
#include "Attributes/Attributes.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/AValue.h"
#include "Utilities/OpalException.h"
#include "Utilities/ParseError.h"
#include <vector>

using namespace Expressions;


// Class StringArray
// ------------------------------------------------------------------------

namespace Attributes {

    StringArray::StringArray(const std::string &name, const std::string &help):
        AttributeHandler(name, help, 0)
    {}


    StringArray::~StringArray()
    {}


    const std::string &StringArray::getType() const {
        static std::string type = "string array";
        return type;
    }


    void StringArray::parse(Attribute &attr, Statement &stat, bool) const {
        Attributes::setStringArray(attr, Expressions::parseStringArray(stat));
    }


    void StringArray::parseComponent
    (Attribute &attr, Statement &statement, bool, int index) const {
        std::vector<std::string> array = Attributes::getStringArray(attr);

        if(AttributeBase *base = &attr.getBase()) {
            array = dynamic_cast<AValue<std::string>*>(base)->evaluate();
        }

        while(int(array.size()) < index) {
            array.push_back(std::string());
        }

        array[index - 1] =
            Expressions::parseStringValue(statement, "String value expected.");
        Attributes::setStringArray(attr, array);
    }

};
