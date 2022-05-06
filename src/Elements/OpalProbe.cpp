//
// Class OpalProbe
//   The Probe element.
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
#include "Elements/OpalProbe.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/ProbeRep.h"
#include "Structure/OpalWake.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"

OpalProbe::OpalProbe():
    OpalElement(SIZE, "PROBE",
                "The \"PROBE\" element defines a Probe."),
    owk_m(nullptr) {

    itsAttr[XSTART] = Attributes::makeReal
                      ("XSTART", " Start of x coordinate [mm]");
    itsAttr[XEND] = Attributes::makeReal
                    ("XEND", " End of x coordinate [mm]");
    itsAttr[YSTART] = Attributes::makeReal
                      ("YSTART", "Start of y coordinate [mm]");
    itsAttr[YEND] = Attributes::makeReal
                    ("YEND", "End of y coordinate [mm]");
    itsAttr[WIDTH] = Attributes::makeReal
                     ("WIDTH", "Width of the probe, not used.");
    itsAttr[STEP] = Attributes::makeReal
                     ("STEP", "Step size of the probe [mm]", 1.0);

    registerOwnership();

    setElement(new ProbeRep("PROBE"));
}


OpalProbe::OpalProbe(const std::string &name, OpalProbe *parent):
    OpalElement(name, parent),
    owk_m(nullptr) {
    setElement(new ProbeRep(name));
}


OpalProbe::~OpalProbe() {
    delete owk_m;
}


OpalProbe *OpalProbe::clone(const std::string &name) {
    return new OpalProbe(name, this);
}


void OpalProbe::update() {
    OpalElement::update();

    ProbeRep *prob =
        dynamic_cast<ProbeRep *>(getElement());

    double xstart = Units::mm2m * Attributes::getReal(itsAttr[XSTART]);
    double xend   = Units::mm2m * Attributes::getReal(itsAttr[XEND]);
    double ystart = Units::mm2m * Attributes::getReal(itsAttr[YSTART]);
    double yend   = Units::mm2m * Attributes::getReal(itsAttr[YEND]);
    double step   = Units::mm2m * Attributes::getReal(itsAttr[STEP]);

    double length = Attributes::getReal(itsAttr[LENGTH]);

    if(itsAttr[WAKEF] && owk_m == nullptr) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + std::string("_wake"));
        owk_m->initWakefunction(*prob);
        prob->setWake(owk_m->wf_m);
    }
    prob->setElementLength(length); // is this needed here?
    prob->setDimensions(xstart, xend, ystart, yend);
    prob->setStep(step);
    prob->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(prob);
}
