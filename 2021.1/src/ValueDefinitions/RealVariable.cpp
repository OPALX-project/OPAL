//
// Class RealVariable
//   The REAL VARIABLE definition.
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
#include "ValueDefinitions/RealVariable.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"

#include <iostream>


RealVariable::RealVariable():
    ValueDefinition(1, "REAL_VARIABLE",
                    "The \"REAL VARIABLE\" statement defines a global "
                    "real variable:\n"
                    "\tREAL VARIABLE <name>=<real-expression>;\n") {
    itsAttr[0] = Attributes::makeReal("VALUE", "The variable value", 0.0);

    registerOwnership(AttributeHandler::STATEMENT);

    // Construct the P0 variable.
    RealVariable *p0 = new RealVariable("P0", this, 1.0);
    OpalData::getInstance()->create(p0);
    p0->setDirty(true);
    OpalData::getInstance()->setP0(p0);
}


RealVariable::RealVariable(const std::string &name, RealVariable *parent,
                           double value):
    ValueDefinition(name, parent) {
    Attributes::setReal(itsAttr[0], value);
}


RealVariable::RealVariable(const std::string &name, RealVariable *parent):
    ValueDefinition(name, parent)
{}


RealVariable::~RealVariable()
{}


bool RealVariable::canReplaceBy(Object *object) {
    // Replace only by another variable.
    return (dynamic_cast<RealVariable *>(object) != 0);
}


RealVariable *RealVariable::clone(const std::string &name) {
    return new RealVariable(name, this);
}


double RealVariable::getReal() const {
    return Attributes::getReal(itsAttr[0]);
}


void RealVariable::print(std::ostream &os) const {
    os << "REAL " << getOpalName()
       << (itsAttr[0].isExpression() ? ":=" : "=") << itsAttr[0] << ';';
    os << std::endl;
}

void RealVariable::printValue(std::ostream &os) const {
    os << itsAttr[0];
}
