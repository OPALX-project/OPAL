//
// Class OpalSlit
//   The SLIT element.
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
#include "Elements/OpalSlit.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/FlexibleCollimatorRep.h"
#include "Structure/ParticleMatterInteraction.h"

#include <boost/regex.hpp>
#include <cstdlib>

OpalSlit::OpalSlit():
    OpalElement(SIZE, "SLIT",
                "The \"SLIT\" element defines a slit."),
    parmatint_m(nullptr) {
    itsAttr[XSIZE] = Attributes::makeReal
                     ("XSIZE", "Horizontal half-aperture in m");
    itsAttr[YSIZE] = Attributes::makeReal
                     ("YSIZE", "Vertical half-aperture in m");

    registerOwnership();

    setElement(new FlexibleCollimatorRep("SLIT"));
}


OpalSlit::OpalSlit(const std::string& name, OpalSlit* parent):
    OpalElement(name, parent),
    parmatint_m(nullptr) {
    setElement(new FlexibleCollimatorRep(name));
}


OpalSlit::~OpalSlit() {
    if (parmatint_m)
        delete parmatint_m;
}


OpalSlit* OpalSlit::clone(const std::string& name) {
    return new OpalSlit(name, this);
}


void OpalSlit::update() {
    OpalElement::update();

    FlexibleCollimatorRep* coll =
        dynamic_cast<FlexibleCollimatorRep*>(getElement());

    double length = Attributes::getReal(itsAttr[LENGTH]);
    coll->setElementLength(length);

    if (getOpalName() != "SLIT") {
        double width = 2 * Attributes::getReal(itsAttr[XSIZE]);
        double height = 2 * Attributes::getReal(itsAttr[YSIZE]);
        std::stringstream description;
        description << "rectangle(" << width << "," << height << ")";
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

    std::vector<double> apert = {Attributes::getReal(itsAttr[XSIZE]),
                                 Attributes::getReal(itsAttr[YSIZE]),
                                 1.0};
    coll->setAperture(ApertureType::RECTANGULAR, apert );

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(coll);
}
