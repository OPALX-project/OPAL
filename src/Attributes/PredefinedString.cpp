//
// Class PredefinedString
//   This class is used to parse attributes of type string that should
//   be contained in set of predefined strings.
//
// Copyright (c) 2021, Christof Metzger-Kraus, Open Sourcerer
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
#include "Attributes/Attributes.h"
#include "Attributes/PredefinedString.h"
#include "AbstractObjects/Expressions.h"
#include "Utilities/ParseError.h"
#include "Utilities/Util.h"

using namespace Expressions;


namespace Attributes {

    PredefinedString::PredefinedString(const std::string &name,
                                       const std::string &help,
                                       const std::initializer_list<std::string>& predefinedStrings,
                                       const std::string& defaultValue):
        AttributeHandler(name, help, 0)
    {
        for (const std::string& value: predefinedStrings) {
            predefinedStrings_m.insert(Util::toUpper(value));
        }
        setPredefinedValues(predefinedStrings_m, defaultValue);
    }


    PredefinedString::~PredefinedString()
    {}


    const std::string &PredefinedString::getType() const {
        static const std::string type("string");
        return type;
    }

    void PredefinedString::parse(Attribute &attr, Statement &stat, bool) const {
        std::string value = Util::toUpper(parseStringValue(stat, "String value expected"));
        if (predefinedStrings_m.count(value) == 0) {
            throw ParseError("PredefinedString::parse", "Unsupported value '" + value + "'");
        }
        Attributes::setPredefinedString(attr, value);
    }

};