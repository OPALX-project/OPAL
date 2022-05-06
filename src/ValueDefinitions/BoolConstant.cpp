//
// Class BoolConstant
//   The BOOL CONSTANT definition.
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
#include "ValueDefinitions/BoolConstant.h"

#include "Attributes/Attributes.h"
#include "Utilities/Options.h"

#include <iostream>


BoolConstant::BoolConstant():
    ValueDefinition(1, "BOOL_CONSTANT",
                    "The \"BOOL CONSTANT\" statement defines a global "
                    "logical constant:\n"
                    "\tBOOL CONSTANT <name> = <Bool-expression>;\n") {
    itsAttr[0] = Attributes::makeBool("VALUE", "The constant value");

    registerOwnership(AttributeHandler::STATEMENT);
}


BoolConstant::BoolConstant(const std::string &name, BoolConstant *parent):
    ValueDefinition(name, parent)
{}


BoolConstant::~BoolConstant()
{}


bool BoolConstant::canReplaceBy(Object *) {
    return false;
}


BoolConstant *BoolConstant::clone(const std::string &name) {
    return new BoolConstant(name, this);
}


bool BoolConstant::getBool() const {
    return Attributes::getBool(itsAttr[0]);
}


void BoolConstant::print(std::ostream &os) const {
    os << "BOOL CONST " << getOpalName() << '=' << itsAttr[0] << ';';
    os << std::endl;
}

void BoolConstant::printValue(std::ostream &os) const {
    os << itsAttr[0];
}
