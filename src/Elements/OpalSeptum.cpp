//
// Class OpalSeptum
//   The Septum element.
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
#include "Elements/OpalSeptum.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/SeptumRep.h"
#include "Structure/OpalWake.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"

OpalSeptum::OpalSeptum():
    OpalElement(SIZE, "SEPTUM",
                "The \"SEPTUM\" element defines a Septum."),
    owk_m(nullptr) {

    itsAttr[XSTART] = Attributes::makeReal
                      ("XSTART", " Start of x coordinate");
    itsAttr[XEND] = Attributes::makeReal
                    ("XEND", " End of x coordinate");
    itsAttr[YSTART] = Attributes::makeReal
                      ("YSTART", "Start of y coordinate");
    itsAttr[YEND] = Attributes::makeReal
                    ("YEND", "End of y coordinate");
    itsAttr[WIDTH] = Attributes::makeReal
                     ("WIDTH", "Width of the septum");

    registerOwnership();

    setElement(new SeptumRep("SEPTUM"));
}


OpalSeptum::OpalSeptum(const std::string &name, OpalSeptum *parent):
    OpalElement(name, parent),
    owk_m(nullptr) {
    setElement(new SeptumRep(name));
}


OpalSeptum::~OpalSeptum() {
    if(owk_m)
        delete owk_m;
}


OpalSeptum *OpalSeptum::clone(const std::string &name) {
    return new OpalSeptum(name, this);
}


void OpalSeptum::update() {
    OpalElement::update();

    SeptumRep *sept = dynamic_cast<SeptumRep *>(getElement());

    double xstart = Units::mm2m * Attributes::getReal(itsAttr[XSTART]);
    double xend   = Units::mm2m * Attributes::getReal(itsAttr[XEND]);
    double ystart = Units::mm2m * Attributes::getReal(itsAttr[YSTART]);
    double yend   = Units::mm2m * Attributes::getReal(itsAttr[YEND]);
    double width  = Units::mm2m * Attributes::getReal(itsAttr[WIDTH]);

    double length = Attributes::getReal(itsAttr[LENGTH]);

    if(itsAttr[WAKEF] && owk_m == nullptr) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + std::string("_wake"));
        owk_m->initWakefunction(*sept);
        sept->setWake(owk_m->wf_m);
    }
    sept->setElementLength(length);
    sept->setDimensions(xstart, xend, ystart, yend);
    sept->setWidth(width);
    sept->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(sept);
}
