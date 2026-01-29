//
// Class OpalSextupole
//   The SEXTUPOLE element.
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
#include "Elements/OpalSextupole.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MultipoleRep.h"
#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include <cmath>
#include <iostream>
#include <sstream>


OpalSextupole::OpalSextupole():
    OpalElement(SIZE, "SEXTUPOLE",
                "The \"SEXTUPOLE\" element defines a Sextupole.") {
    itsAttr[K2] = Attributes::makeReal
                  ("K2", "Normalised upright sextupole coefficient in m^(-3)");
    itsAttr[DK2] = Attributes::makeReal
                  ("DK2", "Normalised upright sextupole coefficient error in m^(-3)");
    itsAttr[K2S] = Attributes::makeReal
                   ("K2S", "Normalised skew sextupole coefficient in m^(-3)");
    itsAttr[DK2S] = Attributes::makeReal
                   ("DK2S", "Normalised skew sextupole coefficient error in m^(-3)");

    registerOwnership();

    setElement(new MultipoleRep("SEXTUPOLE"));
}


OpalSextupole::OpalSextupole(const std::string &name, OpalSextupole *parent):
    OpalElement(name, parent) {
    setElement(new MultipoleRep(name));
}


OpalSextupole::~OpalSextupole()
{}


OpalSextupole *OpalSextupole::clone(const std::string &name) {
    return new OpalSextupole(name, this);
}



void OpalSextupole::print(std::ostream &os) const {
    OpalElement::print(os);
}


void OpalSextupole::update() {
    OpalElement::update();

    MultipoleRep *sext =
        dynamic_cast<MultipoleRep *>(getElement());
    sext->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    double factor = OpalData::getInstance()->getP0() / (Physics::c * 2.0);
    BMultipoleField field;
    field.setNormalComponent(3, factor * Attributes::getReal(itsAttr[K2]));
    field.setSkewComponent(3, factor * Attributes::getReal(itsAttr[K2S]));
    sext->setField(field);
    sext->setNormalComponent(3, Attributes::getReal(itsAttr[K2]), Attributes::getReal(itsAttr[DK2]));
    sext->setSkewComponent(3, Attributes::getReal(itsAttr[K2S]), Attributes::getReal(itsAttr[DK2S]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(sext);
}