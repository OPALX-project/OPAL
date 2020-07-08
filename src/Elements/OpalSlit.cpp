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
    parmatint_m(NULL) {
    itsAttr[XSIZE] = Attributes::makeReal
                     ("XSIZE", "Horizontal half-aperture in m");
    itsAttr[YSIZE] = Attributes::makeReal
                     ("YSIZE", "Vertical half-aperture in m");
    itsAttr[OUTFN] = Attributes::makeString
                     ("OUTFN", "Monitor output filename");


    registerStringAttribute("OUTFN");
    registerRealAttribute("XSIZE");
    registerRealAttribute("YSIZE");

    registerOwnership();

    setElement(new FlexibleCollimatorRep("SLIT"));
}


OpalSlit::OpalSlit(const std::string &name, OpalSlit *parent):
    OpalElement(name, parent),
    parmatint_m(NULL) {
    setElement(new FlexibleCollimatorRep(name));
}


OpalSlit::~OpalSlit() {
    if(parmatint_m)
        delete parmatint_m;
}


OpalSlit *OpalSlit::clone(const std::string &name) {
    return new OpalSlit(name, this);
}


void OpalSlit::fillRegisteredAttributes(const ElementBase &base) {
    OpalElement::fillRegisteredAttributes(base);

    const FlexibleCollimatorRep *coll =
        dynamic_cast<const FlexibleCollimatorRep *>(&base);

    std::string Double("(-?[0-9]+\\.?[0-9]*([Ee][+-]?[0-9]+)?)");
    std::string desc = coll->getDescription();

    boost::regex parser("rectangle\\(" + Double + "," + Double + "\\)");
    boost::smatch what;
    if (!boost::regex_match(desc, what, parser)) return;

    double width = atof(std::string(what[1]).c_str());
    double height = atof(std::string(what[3]).c_str());

    attributeRegistry["XSIZE"]->setReal(0.5 * width);
    attributeRegistry["YSIZE"]->setReal(0.5 * height);
}


void OpalSlit::update() {
    OpalElement::update();

    FlexibleCollimatorRep *coll =
        dynamic_cast<FlexibleCollimatorRep *>(getElement());
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

    if(itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == NULL) {
        parmatint_m = (ParticleMatterInteraction::find(Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION])))->clone(getOpalName() + std::string("_parmatint"));
        parmatint_m->initParticleMatterInteractionHandler(*coll);
        coll->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    std::vector<double> apert = {Attributes::getReal(itsAttr[XSIZE]),
                                 Attributes::getReal(itsAttr[YSIZE]),
                                 1.0};
    coll->setAperture(ElementBase::RECTANGULAR, apert );

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(coll);
}