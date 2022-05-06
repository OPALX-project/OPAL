//
// Class OpalHKicker
//   The HKICKER element.
//   Note the sign convention:  A positive kick bend particles to positive x.
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
#include "Elements/OpalHKicker.h"

#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/XCorrectorRep.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"

OpalHKicker::OpalHKicker():
    OpalElement(SIZE, "HKICKER",
                "The \"HKICKER\" element defines a closed orbit corrector "
                "acting on the horizontal plane.") {
    itsAttr[KICK] = Attributes::makeReal
                    ("KICK", "Horizontal deflection in rad");
    itsAttr[DESIGNENERGY] = Attributes::makeReal
                           ("DESIGNENERGY", "the mean energy of the particles");
    itsAttr[K0] = Attributes::makeReal
                  ("K0", "Normal dipole field in T");

    registerOwnership();

    setElement(new XCorrectorRep("HKICKER"));
}


OpalHKicker::OpalHKicker(const std::string &name, OpalHKicker *parent):
    OpalElement(name, parent) {
    setElement(new XCorrectorRep(name));
}


OpalHKicker::~OpalHKicker()
{}


OpalHKicker *OpalHKicker::clone(const std::string &name) {
    return new OpalHKicker(name, this);
}


void OpalHKicker::update() {
    OpalElement::update();

    XCorrectorRep *corr =
        dynamic_cast<XCorrectorRep *>(getElement());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double factor = OpalData::getInstance()->getP0() / Physics::c;
    double kick = Attributes::getReal(itsAttr[KICK]);

    corr->setElementLength(length);
    corr->setBy(- kick * factor);

    corr->setKickX(kick);
    if(itsAttr[DESIGNENERGY]) {
        double kineticEnergy = Attributes::getReal(itsAttr[DESIGNENERGY]) * Units::MeV2eV;
        corr->setDesignEnergy(kineticEnergy, false);
    }

    if (itsAttr[K0]) {
        corr->setKickField(Vector_t(0, Attributes::getReal(itsAttr[K0]), 0));
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(corr);
}