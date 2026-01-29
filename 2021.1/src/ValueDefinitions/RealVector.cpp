//
// Class RealVector
//   The REAL VECTOR definition.
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
#include "ValueDefinitions/RealVector.h"

#include "Attributes/Attributes.h"
#include "Utilities/Options.h"

#include <iostream>
#include <vector>


RealVector::RealVector():
    ValueDefinition(1, "REAL_VECTOR",
                    "The \"REAL VECTOR\" statement defines a global "
                    "real vector.\n"
                    "\tREAL VECTOR<name>=<real-vector-expression>;\n") {
    itsAttr[0] = Attributes::makeRealArray("VALUE", "The vector value");

    registerOwnership(AttributeHandler::STATEMENT);
}


RealVector::RealVector(const std::string &name, RealVector *parent):
    ValueDefinition(name, parent)
{}


RealVector::~RealVector()
{}


bool RealVector::canReplaceBy(Object *object) {
    // Replace only by another vector.
    return (dynamic_cast<RealVector *>(object) != 0);
}


RealVector *RealVector::clone(const std::string &name) {
    return new RealVector(name, this);
}


void RealVector::print(std::ostream &os) const {
    // WARNING: Cannot print in OPAL-8 format.
    os << "REAL VECTOR " << getOpalName() << ":="
       << itsAttr[0] << ';' << std::endl;
}

void RealVector::printValue(std::ostream &os) const {
    os << itsAttr[0];
}

double RealVector::getRealComponent(int index) const {
    std::vector<double> array = Attributes::getRealArray(itsAttr[0]);
    return array[index-1];
}
