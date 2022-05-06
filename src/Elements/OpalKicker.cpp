//
// Class OpalKicker
//   The KICKER element.
//   Note the sign convention:  Positive kicks bend particles to positive x or
//   y respectively.
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
#include "Elements/OpalKicker.h"

#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/CorrectorRep.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"


OpalKicker::OpalKicker():
    OpalElement(SIZE, "KICKER",
                "The \"KICKER\" element defines a closed orbit corrector "
                "acting on both planes.") {
    itsAttr[HKICK] = Attributes::makeReal
                     ("HKICK", "Horizontal deflection in rad");
    itsAttr[VKICK] = Attributes::makeReal
                     ("VKICK", "Vertical deflection in rad");
    itsAttr[DESIGNENERGY] = Attributes::makeReal
                            ("DESIGNENERGY", "the mean energy of the particles");
    itsAttr[K0] = Attributes::makeReal
                  ("K0", "Normal dipole field in T");
    itsAttr[K0S] = Attributes::makeReal
                  ("K0S", "Skew dipole field in T");

    registerOwnership();

    setElement(new CorrectorRep("KICKER"));
}


OpalKicker::OpalKicker(const std::string &name, OpalKicker *parent):
    OpalElement(name, parent) {
    setElement(new CorrectorRep(name));
}


OpalKicker::~OpalKicker()
{}


OpalKicker *OpalKicker::clone(const std::string &name) {
    return new OpalKicker(name, this);
}


void OpalKicker::update() {
    OpalElement::update();

    CorrectorRep *corr =
        dynamic_cast<CorrectorRep *>(getElement());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double factor = OpalData::getInstance()->getP0() / Physics::c;
    double hKick = Attributes::getReal(itsAttr[HKICK]);
    double vKick = Attributes::getReal(itsAttr[VKICK]);

    corr->setElementLength(length);
    corr->setBy(- hKick * factor);
    corr->setBx(vKick * factor);

    corr->setKickX(hKick);
    corr->setKickY(vKick);
    if(itsAttr[DESIGNENERGY]) {
        double kineticEnergy = Attributes::getReal(itsAttr[DESIGNENERGY]) * Units::MeV2eV;
        corr->setDesignEnergy(kineticEnergy, false);
    }

    double Bx = 0.0, By = 0.0;
    bool fieldSet = false;
    if (itsAttr[K0]) {
        Bx = Attributes::getReal(itsAttr[K0]);
        fieldSet = true;
    }
    if (itsAttr[K0S]) {
        By = Attributes::getReal(itsAttr[K0S]);
        fieldSet = true;
    }

    if (fieldSet) {
        corr->setKickField(Vector_t(Bx, By, 0));
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(corr);
}