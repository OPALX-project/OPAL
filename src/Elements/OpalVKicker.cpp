//
// Class OpalVKicker
//   The VKICKER element.
//   Note the sign convention:  A positive kick bend particles to positive y.
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
#include "Elements/OpalVKicker.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/YCorrectorRep.h"
#include "Physics/Physics.h"


OpalVKicker::OpalVKicker():
    OpalElement(SIZE, "VKICKER",
                "The \"VKICKER\" element defines a closed orbit corrector "
                "acting on the vertical plane.") {
    itsAttr[KICK] = Attributes::makeReal
                    ("KICK", "Vertical deflection in rad");
    itsAttr[DESIGNENERGY] = Attributes::makeReal
                            ("DESIGNENERGY", "the mean energy of the particles");
    itsAttr[K0] = Attributes::makeReal
                  ("K0", "Normal dipole field in T");

    registerRealAttribute("KICK");
    registerRealAttribute("DESIGNENERGY");
    registerRealAttribute("K0");

    registerOwnership();

    setElement(new YCorrectorRep("VKICKER"));
}


OpalVKicker::OpalVKicker(const std::string &name, OpalVKicker *parent):
    OpalElement(name, parent) {
    setElement(new YCorrectorRep(name));
}


OpalVKicker::~OpalVKicker()
{}


OpalVKicker *OpalVKicker::clone(const std::string &name) {
    return new OpalVKicker(name, this);
}


void OpalVKicker::
fillRegisteredAttributes(const ElementBase &base) {
    Inform m("fillRegisteredAttributes ");
    OpalElement::fillRegisteredAttributes(base);

    const YCorrectorRep *corr = dynamic_cast<const YCorrectorRep *>(&base);
    BDipoleField field = corr->getField();

    double scale = Physics::c / OpalData::getInstance()->getP0();
    attributeRegistry["HKICK"]->setReal(- field.getBy() * scale);
    attributeRegistry["VKICK"]->setReal(+ field.getBx() * scale);

    m << "P= " << OpalData::getInstance()->getP0()
      << " Bx= " << field.getBx()
      << " By= " << field.getBy() << endl;

}


void OpalVKicker::update() {
    OpalElement::update();

    YCorrectorRep *corr = dynamic_cast<YCorrectorRep *>(getElement());

    double length = Attributes::getReal(itsAttr[LENGTH]);
    double factor = OpalData::getInstance()->getP0() / Physics::c;
    double kick = Attributes::getReal(itsAttr[KICK]);

    corr->setElementLength(length);
    corr->setBx(kick * factor);
    corr->setKickY(kick);
    if(itsAttr[DESIGNENERGY]) {
        double kineticEnergy = Attributes::getReal(itsAttr[DESIGNENERGY]) * 1e6;
        corr->setDesignEnergy(kineticEnergy, false);
    }

    if (itsAttr[K0]) {
        corr->setKickField(Vector_t(Attributes::getReal(itsAttr[K0]), 0, 0));
    }
    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(corr);
}