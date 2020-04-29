//
// Class UpperCaseStringArray
//   This class is used to parse attributes of type array of strings that should
//   be all upper case, i.e. it returns an array of upper case strings.
//
// Copyright (c) 2020, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Attributes/UpperCaseStringArray.h"
#include "Attributes/Attributes.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/AValue.h"
#include "Utilities/OpalException.h"
#include "Utilities/ParseError.h"
#include "Utilities/Util.h"
#include <vector>

using namespace Expressions;


namespace Attributes {

    UpperCaseStringArray::UpperCaseStringArray(const std::string &name, const std::string &help):
        AttributeHandler(name, help, 0)
    {}


    UpperCaseStringArray::~UpperCaseStringArray()
    {}


    const std::string &UpperCaseStringArray::getType() const {
        static std::string type = "upper case string array";
        return type;
    }


    void UpperCaseStringArray::parse(Attribute &attr, Statement &stat, bool) const {
        std::vector<std::string> array = Expressions::parseStringArray(stat);
        Attributes::setUpperCaseStringArray(attr, array);
    }


    void UpperCaseStringArray::parseComponent
    (Attribute &attr, Statement &statement, bool, int index) const {
        std::vector<std::string> array = Attributes::getStringArray(attr);

        if(AttributeBase *base = &attr.getBase()) {
            array = dynamic_cast<AValue<std::string>*>(base)->evaluate();
        }

        while(int(array.size()) < index) {
            array.push_back(std::string());
        }

        array[index - 1] =
            Expressions::parseString(statement, "String value expected.");
        Attributes::setUpperCaseStringArray(attr, array);
    }

};
