//
// Class Value
//   The class for OPAL VALUE command.
//
// Copyright (c) 2000 - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "BasicActions/Value.h"

#include "Attributes/Attributes.h"
#include "Parser/Statement.h"
#include "Utility/Inform.h"

#include <string>
#include <vector>


extern Inform* gmsg;

Value::Value():
    Action(1, "VALUE",
           "The \"VALUE\" statement prints a list of expressions and "
           "their values.") {
    itsAttr[0] = Attributes::makeRealArray
                 ("VALUE", "The values to be evaluated");

    registerOwnership(AttributeHandler::STATEMENT);
}


Value::Value(const std::string& name, Value* parent):
    Action(name, parent)
{}


Value::~Value()
{}


Value* Value::clone(const std::string& name) {
    return new Value(name, this);
}


void Value::execute() {
    *gmsg << "\nvalue: " << itsAttr[0] << " = {";
    //  std::streamsize old_prec = *gmsg.precision(12);
    const std::vector<double> array = Attributes::getRealArray(itsAttr[0]);
    std::vector<double>::const_iterator i = array.begin();

    while (i != array.end()) {
        *gmsg << *i++;
        if (i == array.end()) break;
        *gmsg << ", ";
    }

    *gmsg << "}\n" << endl;
    //  *gmsg.precision(old_prec);
}

void Value::parse(Statement& statement) {
    // parse, but don't evaluate (for printing mainly)
    parseShortcut(statement, false);
}
