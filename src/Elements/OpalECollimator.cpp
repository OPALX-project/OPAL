//
// Class OpalECollimator
//   The ECOLLIMATOR element.
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
#include "Elements/OpalECollimator.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/FlexibleCollimatorRep.h"
#include "Structure/ParticleMatterInteraction.h"

#include <boost/regex.hpp>
#include <cstdlib>

OpalECollimator::OpalECollimator():
    OpalElement(SIZE, "ECOLLIMATOR",
                "The \"ECOLLIMATOR\" element defines an elliptic collimator."),
    parmatint_m(nullptr) {
    itsAttr[XSIZE] = Attributes::makeReal
                     ("XSIZE", "Horizontal half-aperture in m");
    itsAttr[YSIZE] = Attributes::makeReal
                     ("YSIZE", "Vertical half-aperture in m");

    registerOwnership();

    setElement(new FlexibleCollimatorRep("ECOLLIMATOR"));
}


OpalECollimator::OpalECollimator(const std::string& name, OpalECollimator* parent):
    OpalElement(name, parent),
    parmatint_m(nullptr) {
    setElement(new FlexibleCollimatorRep(name));
}


OpalECollimator::~OpalECollimator() {
    if (parmatint_m)
        delete parmatint_m;
}


OpalECollimator* OpalECollimator::clone(const std::string& name) {
    return new OpalECollimator(name, this);
}


void OpalECollimator::update() {
    OpalElement::update();

    FlexibleCollimatorRep* coll =
        dynamic_cast<FlexibleCollimatorRep*>(getElement());

    double length = Attributes::getReal(itsAttr[LENGTH]);
    coll->setElementLength(length);

    if (getOpalName() != "ECOLLIMATOR") {
        double width  = 2 * Attributes::getReal(itsAttr[XSIZE]);
        double height = 2 * Attributes::getReal(itsAttr[YSIZE]);
        std::stringstream description;
        description << "ellipse(" << width << "," << height << ")";
        coll->setDescription(description.str());
    }

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
