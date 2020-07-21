//
// Class OpalPepperPot
//   The PEPPERPOT element.
//   The class of OPAL elliptic collimators.
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
#include "Elements/OpalPepperPot.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/FlexibleCollimatorRep.h"
#include "Structure/ParticleMatterInteraction.h"


OpalPepperPot::OpalPepperPot():
    OpalElement(SIZE, "PEPPERPOT",
                "The \"PEPPERPOT\" element defines an elliptic collimator."),
    parmatint_m(NULL) {
    itsAttr[XSIZE] = Attributes::makeReal
                     ("XSIZE", "Size in x of the pepperpot in m");
    itsAttr[YSIZE] = Attributes::makeReal
                     ("YSIZE", "Size in y of the pepperpot in m");
    itsAttr[OUTFN] = Attributes::makeString
                     ("OUTFN", "Pepperpot output filename");
    itsAttr[NHOLX] = Attributes::makeReal
                     ("NHOLX", "Number of holes in x");
    itsAttr[NHOLY] = Attributes::makeReal
                     ("NHOLY", "Number of holes in y");
    itsAttr[R] = Attributes::makeReal
                 ("R", "Radios of a holes in m");

    registerOwnership();

    setElement(new FlexibleCollimatorRep("PEPPERPOT"));
}


OpalPepperPot::OpalPepperPot(const std::string &name, OpalPepperPot *parent):
    OpalElement(name, parent),
    parmatint_m(NULL) {
    setElement(new FlexibleCollimatorRep(name));
}


OpalPepperPot::~OpalPepperPot() {
    if(parmatint_m)
        delete parmatint_m;
}


OpalPepperPot *OpalPepperPot::clone(const std::string &name) {
    return new OpalPepperPot(name, this);
}


void OpalPepperPot::update() {
    OpalElement::update();

    FlexibleCollimatorRep *ppo =
        dynamic_cast<FlexibleCollimatorRep *>(getElement());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    ppo->setElementLength(length);
    ppo->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    if (getOpalName() != "PEPPERPOT") {
        double xsize = Attributes::getReal(itsAttr[XSIZE]);
        double ysize = Attributes::getReal(itsAttr[YSIZE]);
        double diameter = 2 * Attributes::getReal(itsAttr[R]);
        int repX = Attributes::getReal(itsAttr[NHOLX]) - 1;
        int repY = Attributes::getReal(itsAttr[NHOLY]) - 1;

        double shiftx = (xsize - diameter) / repX;
        double shifty = (ysize - diameter) / repY;

        std::stringstream description;
        description << "repeat(repeat(translate(ellipse("
                    << diameter << "," << diameter << "),"
                    << -shiftx * 0.5 * repX << "," << -shifty * 0.5 * repY << "),"
                    << repX << "," << shiftx << ",0.0),"
                    << repY << ",0.0," << shifty << ")";

        std::cout << "OpalPepperPot.cpp: " << __LINE__ << "\t"
                  << description.str() << std::endl;
        ppo->setDescription(description.str());
        exit(1);
    }

    if(itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == NULL) {
        parmatint_m = (ParticleMatterInteraction::find(Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION])))->clone(getOpalName() + std::string("_parmatint"));
        parmatint_m->initParticleMatterInteractionHandler(*ppo);
        ppo->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(ppo);
}