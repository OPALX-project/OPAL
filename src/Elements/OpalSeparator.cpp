//
// Class OpalSeparator
//   The ELSEPARATOR element.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Elements/OpalSeparator.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/SeparatorRep.h"


OpalSeparator::OpalSeparator():
    OpalElement(SIZE, "SEPARATOR",
                "The \"SEPARATOR\" element defines an electrostatic separator.") {
    itsAttr[EX] = Attributes::makeReal
                  ("EX", "The horizontal electrostatic field in MV");
    itsAttr[EY] = Attributes::makeReal
                  ("EY", "The vertical electrostatic field in MV");

    registerRealAttribute("EXL");
    registerRealAttribute("EYL");

    registerOwnership();

    setElement(new SeparatorRep("SEPARATOR"));
}


OpalSeparator::OpalSeparator(const std::string &name, OpalSeparator *parent):
    OpalElement(name, parent) {
    setElement(new SeparatorRep(name));
}


OpalSeparator::~OpalSeparator()
{}


OpalSeparator *OpalSeparator::clone(const std::string &name) {
    return new OpalSeparator(name, this);
}


void OpalSeparator::
fillRegisteredAttributes(const ElementBase &base) {
    OpalElement::fillRegisteredAttributes(base);

    const SeparatorRep *sep = dynamic_cast<const SeparatorRep *>(&base);
    double length = sep->getElementLength();
    attributeRegistry["EXL"]->setReal(length * sep->getEx());
    attributeRegistry["EYL"]->setReal(length * sep->getEy());
}


void OpalSeparator::update() {
    OpalElement::update();

    SeparatorRep *sep =
        dynamic_cast<SeparatorRep *>(getElement());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double Ex     = Attributes::getReal(itsAttr[EX]) * 1.0e6;
    double Ey     = Attributes::getReal(itsAttr[EY]) * 1.0e6;
    sep->setElementLength(length);
    sep->setEx(Ex);
    sep->setEy(Ey);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(sep);
}