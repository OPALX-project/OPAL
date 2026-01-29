//
// Class OpalOctupole
//   The OCTUPOLE element.
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
#include "Elements/OpalOctupole.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MultipoleRep.h"
#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include <cmath>
#include <iostream>
#include <sstream>


OpalOctupole::OpalOctupole():
    OpalElement(SIZE, "OCTUPOLE",
                "The \"OCTUPOLE\" element defines a Octupole.") {
    itsAttr[K3] = Attributes::makeReal
                  ("K3", "Normalised upright octupole coefficient in m^(-4)");
    itsAttr[DK3] = Attributes::makeReal
                  ("DK3", "Normalised upright octupole coefficient error in m^(-4)");
    itsAttr[K3S] = Attributes::makeReal
                   ("K3S", "Normalised skew octupole coefficient in m^(-4)");
    itsAttr[DK3S] = Attributes::makeReal
                   ("DK3S", "Normalised skew octupole coefficient error in m^(-4)");

    registerOwnership();

    setElement(new MultipoleRep("OCTUPOLE"));
}


OpalOctupole::OpalOctupole(const std::string &name, OpalOctupole *parent):
    OpalElement(name, parent) {
    setElement(new MultipoleRep(name));
}


OpalOctupole::~OpalOctupole()
{}


OpalOctupole *OpalOctupole::clone(const std::string &name) {
    return new OpalOctupole(name, this);
}


void OpalOctupole::print(std::ostream &os) const {
    OpalElement::print(os);
}


void OpalOctupole::update() {
    OpalElement::update();

    MultipoleRep *oct =
        dynamic_cast<MultipoleRep *>(getElement());
    oct->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    double factor = OpalData::getInstance()->getP0() / (Physics::c * 6.0);
    BMultipoleField field;
    field.setNormalComponent(4, factor * Attributes::getReal(itsAttr[K3]));
    field.setSkewComponent(4, factor * Attributes::getReal(itsAttr[K3S]));
    oct->setField(field);

    oct->setNormalComponent(4, Attributes::getReal(itsAttr[K3]), Attributes::getReal(itsAttr[DK3]));
    oct->setSkewComponent(4, Attributes::getReal(itsAttr[K3S]), Attributes::getReal(itsAttr[DK3S]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(oct);
}