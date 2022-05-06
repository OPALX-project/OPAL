//
// Class OpalCCollimator
//   The CCOLLIMATOR element.
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
#include "Elements/OpalCCollimator.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/CCollimatorRep.h"
#include "Structure/ParticleMatterInteraction.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"

OpalCCollimator::OpalCCollimator():
    OpalElement(SIZE, "CCOLLIMATOR",
                "The \"CCOLLIMATOR\" element defines a rectangular-shape cyclotron collimator"),
    parmatint_m(nullptr) {
    itsAttr[XSTART] = Attributes::makeReal
                      ("XSTART", " Start of x coordinate [mm]");
    itsAttr[XEND]   = Attributes::makeReal
                      ("XEND", " End of x coordinate, [mm]");
    itsAttr[YSTART] = Attributes::makeReal
                      ("YSTART", "Start of y coordinate, [mm]");
    itsAttr[YEND]   = Attributes::makeReal
                      ("YEND", "End of y coordinate, [mm]");
    itsAttr[ZSTART] = Attributes::makeReal
                      ("ZSTART", "Start of vertical coordinate, [mm], default value: -100",-100.0);
    itsAttr[ZEND]   = Attributes::makeReal
                      ("ZEND", "End of vertical coordinate, [mm], default value: 100", 100.0);
    itsAttr[WIDTH]  = Attributes::makeReal
                      ("WIDTH", "Width of the collimator [mm]");

    registerOwnership();

    setElement(new CCollimatorRep("CCOLLIMATOR"));
}


OpalCCollimator::OpalCCollimator(const std::string& name, OpalCCollimator* parent):
    OpalElement(name, parent),
    parmatint_m(nullptr) {
    setElement(new CCollimatorRep(name));
}


OpalCCollimator::~OpalCCollimator() {
    if (parmatint_m)
        delete parmatint_m;
}


OpalCCollimator* OpalCCollimator::clone(const std::string& name) {
    return new OpalCCollimator(name, this);
}


void OpalCCollimator::update() {
    OpalElement::update();

    CCollimatorRep* coll =
        dynamic_cast<CCollimatorRep*>(getElement());

    double xstart = Units::mm2m * Attributes::getReal(itsAttr[XSTART]);
    double xend   = Units::mm2m * Attributes::getReal(itsAttr[XEND]);
    double ystart = Units::mm2m * Attributes::getReal(itsAttr[YSTART]);
    double yend   = Units::mm2m * Attributes::getReal(itsAttr[YEND]);
    double zstart = Units::mm2m * Attributes::getReal(itsAttr[ZSTART]);
    double zend   = Units::mm2m * Attributes::getReal(itsAttr[ZEND]);
    double width  = Units::mm2m * Attributes::getReal(itsAttr[WIDTH]);

    double length = Attributes::getReal(itsAttr[LENGTH]);

    coll->setElementLength(length);
    coll->setDimensions(xstart, xend, ystart, yend, zstart, zend, width);
    coll->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    if (itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == nullptr) {
        const std::string matterDescriptor = Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION]);
        ParticleMatterInteraction* orig = ParticleMatterInteraction::find(matterDescriptor);
        parmatint_m = orig->clone(matterDescriptor);
        parmatint_m->initParticleMatterInteractionHandler(*coll);
        coll->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(coll);
}
