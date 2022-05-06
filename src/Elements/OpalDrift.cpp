//
// Class OpalDrift
//   The class of OPAL drift spaces.
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
#include "Elements/OpalDrift.h"
#include "Structure/BoundaryGeometry.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/DriftRep.h"
#include "Structure/OpalWake.h"
#include "Structure/ParticleMatterInteraction.h"


OpalDrift::OpalDrift():
    OpalElement(SIZE, "DRIFT",
                "The \"DRIFT\" element defines a drift space."),
    owk_m(nullptr),
    parmatint_m(nullptr),
    obgeo_m(nullptr) {
    // CKR: the following 3 lines are redundant: OpalElement does this already!
    //      they prevent drift from working properly
    //
    //     itsAttr[LENGTH] = Attributes::makeReal
    //         ("LENGTH", "Drift length");

    itsAttr[GEOMETRY] = Attributes::makeString
                        ("GEOMETRY", "BoundaryGeometry for Drifts");

    itsAttr[NSLICES]  = Attributes::makeReal
                        ("NSLICES", "The number of slices/ steps for this element in Map Tracking", 1);

    registerOwnership();

    setElement(new DriftRep("DRIFT"));
}


OpalDrift::OpalDrift(const std::string& name, OpalDrift* parent):
    OpalElement(name, parent),
    owk_m(nullptr),
    parmatint_m(nullptr),
    obgeo_m(nullptr) {
    setElement(new DriftRep(name));
}


OpalDrift::~OpalDrift() {
    if (owk_m)
        delete owk_m;
    if (parmatint_m)
        delete parmatint_m;
    if (obgeo_m)
        delete obgeo_m;
}


OpalDrift* OpalDrift::clone(const std::string& name) {
    return new OpalDrift(name, this);
}


bool OpalDrift::isDrift() const {
    return true;
}


void OpalDrift::update() {
    OpalElement::update();

    DriftRep* drf = static_cast<DriftRep*>(getElement());

    drf->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    drf->setNSlices(Attributes::getReal(itsAttr[NSLICES]));
    
    if (itsAttr[WAKEF] && owk_m == nullptr) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + std::string("_wake"));
        owk_m->initWakefunction(*drf);
        drf->setWake(owk_m->wf_m);
    }

    if (itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == nullptr) {
        const std::string matterDescriptor = Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION]);
        ParticleMatterInteraction* orig = ParticleMatterInteraction::find(matterDescriptor);
        parmatint_m = orig->clone(matterDescriptor);
        parmatint_m->initParticleMatterInteractionHandler(*drf);
        drf->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    if (itsAttr[GEOMETRY] && obgeo_m == nullptr) {
        obgeo_m = (BoundaryGeometry::find(Attributes::getString(itsAttr[GEOMETRY])))->clone(getOpalName() + std::string("_geometry"));
        if (obgeo_m) {
            drf->setBoundaryGeometry(obgeo_m);
        }
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(drf);
}
