//
// Class OpalStripper
//   The Stripper element
//
// Copyright (c) 2011, Jianjun Yang, Paul Scherrer Institut, Villigen PSI, Switzerland
// Copyright (c) 2014, 2017-2018, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects"
// and the paper
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects:
// Model, implementation, and application"
// (https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.13.064201)
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
#include "Elements/OpalStripper.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/StripperRep.h"
#include "Physics/Physics.h"


OpalStripper::OpalStripper():
    OpalElement(SIZE, "STRIPPER",
                "The \"STRIPPER\" element defines a Stripper.") {

    itsAttr[XSTART] = Attributes::makeReal
                      ("XSTART", " Start of x coordinate [mm]");
    itsAttr[XEND] = Attributes::makeReal
                    ("XEND", " End of x coordinate, [mm]");
    itsAttr[YSTART] = Attributes::makeReal
                      ("YSTART", "Start of y coordinate, [mm]");
    itsAttr[YEND] = Attributes::makeReal
                    ("YEND", "End of y coordinate, [mm]");
    itsAttr[WIDTH] = Attributes::makeReal
                     ("WIDTH", "Width of the stripper [mm], NOT used");
    itsAttr[OPCHARGE] = Attributes::makeReal
                     ("OPCHARGE", "Charge number of the outcome particle");
    itsAttr[OPMASS] = Attributes::makeReal
                     ("OPMASS", "Mass of the outcome particle [GeV/c^2]");
    itsAttr[OPYIELD] = Attributes::makeReal
                     ("OPYIELD", "Yield (Particle number of the outcome particle) per income particle");
    itsAttr[STOP] = Attributes::makeBool
      ("STOP", "Option Whether stop tracking at the stripper. Default: true", true);
    itsAttr[OUTFN] = Attributes::makeString
                     ("OUTFN", "Output filename");

    registerRealAttribute("XSTART");
    registerRealAttribute("XEND");
    registerRealAttribute("YSTART");
    registerRealAttribute("YEND");
    registerRealAttribute("WIDTH");
    registerRealAttribute("OPCHARGE");
    registerRealAttribute("OPMASS");
    registerRealAttribute("OPYIELD");
    registerStringAttribute("OUTFN");

    registerOwnership();

    setElement((new StripperRep("STRIPPER"))->makeAlignWrapper());
}


OpalStripper::OpalStripper(const std::string &name, OpalStripper *parent):
    OpalElement(name, parent) {
    setElement((new StripperRep(name))->makeAlignWrapper());
}


OpalStripper::~OpalStripper()
{}


OpalStripper *OpalStripper::clone(const std::string &name) {
    return new OpalStripper(name, this);
}


void OpalStripper::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

}


void OpalStripper::update() {
    OpalElement::update();

    StripperRep *strp = dynamic_cast<StripperRep *>(getElement()->removeWrappers());

    const double mm2m = 0.001;
    double xstart   = mm2m * Attributes::getReal(itsAttr[XSTART]);
    double xend     = mm2m * Attributes::getReal(itsAttr[XEND]);
    double ystart   = mm2m * Attributes::getReal(itsAttr[YSTART]);
    double yend     = mm2m * Attributes::getReal(itsAttr[YEND]);

    double length   = Attributes::getReal(itsAttr[LENGTH]);
    double opcharge = Attributes::getReal(itsAttr[OPCHARGE]);
    double opmass   = Attributes::getReal(itsAttr[OPMASS]);
    double opyield  = Attributes::getReal(itsAttr[OPYIELD]);
    bool   stop     = Attributes::getBool(itsAttr[STOP]);

    strp->setElementLength(length);
    strp->setDimensions(xstart, xend, ystart, yend);
    strp->setOPCharge(opcharge);
    strp->setOPMass(opmass);
    strp->setOPYield(opyield);
    strp->setStop(stop);
    strp->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(strp);
}