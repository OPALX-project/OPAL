//
// Class OpalQuadrupole
//   The QUADRUPOLE element.
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
#include "Elements/OpalQuadrupole.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MultipoleRep.h"
#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include "Structure/ParticleMatterInteraction.h"

#include <cmath>
#include <iostream>
#include <sstream>


OpalQuadrupole::OpalQuadrupole():
    OpalElement(SIZE, "QUADRUPOLE",
                "The \"QUADRUPOLE\" element defines a Quadrupole."),
    parmatint_m(nullptr) {
    itsAttr[K1]    = Attributes::makeReal
                     ("K1", "Normalised upright quadrupole coefficient in m^(-2)");
    itsAttr[DK1]   = Attributes::makeReal
                     ("DK1", "Normalised upright quadrupole coefficient error in m^(-2)");
    itsAttr[K1S]   = Attributes::makeReal
                     ("K1S", "Normalised skew quadrupole coefficient in m^(-2)");
    itsAttr[DK1S]  = Attributes::makeReal
                     ("DK1S", "Normalised skew quadrupole coefficient error in m^(-2)");
    itsAttr[NSLICES] = Attributes::makeReal
                       ("NSLICES", "The number of slices/ steps for this element in Map Tracking", 1);

    registerOwnership();

    setElement((new MultipoleRep("QUADRUPOLE")));
}


OpalQuadrupole::OpalQuadrupole(const std::string& name, OpalQuadrupole* parent):
    OpalElement(name, parent),
    parmatint_m(nullptr) {
    setElement((new MultipoleRep(name)));
}


OpalQuadrupole::~OpalQuadrupole() {
    delete parmatint_m;
}


OpalQuadrupole* OpalQuadrupole::clone(const std::string& name) {
    return new OpalQuadrupole(name, this);
}


void OpalQuadrupole::print(std::ostream& os) const {
    OpalElement::print(os);
}


void OpalQuadrupole::update() {
    OpalElement::update();

    MultipoleRep* quad =
        dynamic_cast<MultipoleRep*>(getElement());

    quad->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    double factor = OpalData::getInstance()->getP0() / Physics::c;

    BMultipoleField field;
    field.setNormalComponent(2, factor * Attributes::getReal(itsAttr[K1]));   // this is for the maps
    field.setSkewComponent(2, factor * Attributes::getReal(itsAttr[K1S]));    // this is for the maps
    quad->setField(field);
    quad->setNormalComponent(2, Attributes::getReal(itsAttr[K1]), Attributes::getReal(itsAttr[DK1]));
    quad->setSkewComponent(2, Attributes::getReal(itsAttr[K1S]), Attributes::getReal(itsAttr[DK1S]));
    quad->setNSlices(Attributes::getReal(itsAttr[NSLICES]));

    if (itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == nullptr) {
        const std::string matterDescriptor = Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION]);
        ParticleMatterInteraction* orig = ParticleMatterInteraction::find(matterDescriptor);
        parmatint_m = orig->clone(matterDescriptor);
        parmatint_m->initParticleMatterInteractionHandler(*quad);
        quad->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(quad);
}
