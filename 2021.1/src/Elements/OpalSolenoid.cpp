//
// Class OpalSolenoid
//   The SOLENOID element.
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
#include "Elements/OpalSolenoid.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/SolenoidRep.h"
#include "Physics/Physics.h"


OpalSolenoid::OpalSolenoid():
    OpalElement(SIZE, "SOLENOID",
                "The \"SOLENOID\" element defines a Solenoid.") {
    itsAttr[KS] = Attributes::makeReal
                  ("KS", "Normalised solenoid strength in m**(-1)");
    itsAttr[DKS] = Attributes::makeReal
                  ("DKS", "Normalised solenoid strength error in m**(-1)");
    itsAttr[FMAPFN] = Attributes::makeString
                      ("FMAPFN", "Solenoid field map filename ");
    itsAttr[FAST] = Attributes::makeBool
                    ("FAST", "Faster but less accurate", true);

    registerOwnership();

    setElement(new SolenoidRep("SOLENOID"));
}


OpalSolenoid::OpalSolenoid(const std::string &name, OpalSolenoid *parent):
    OpalElement(name, parent) {
    setElement(new SolenoidRep(name));
}


OpalSolenoid::~OpalSolenoid()
{}


OpalSolenoid *OpalSolenoid::clone(const std::string &name) {
    return new OpalSolenoid(name, this);
}


void OpalSolenoid::update() {
    OpalElement::update();

    SolenoidRep *sol =
        dynamic_cast<SolenoidRep *>(getElement());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double Bz = Attributes::getReal(itsAttr[KS]) * OpalData::getInstance()->getP0() / Physics::c;
    bool fast = Attributes::getBool(itsAttr[FAST]);

    sol->setElementLength(length);
    sol->setFieldMapFN(Attributes::getString(itsAttr[FMAPFN]));
    sol->setFast(fast);
    sol->setBz(Bz);
    sol->setKS(Attributes::getReal(itsAttr[KS]));
    sol->setDKS(Attributes::getReal(itsAttr[DKS]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(sol);
}