//
// Class OpalMultipoleT
//   The Opal MultipoleT element.
//
// Copyright (c) 2017 - 2023, Titus Dascalu
//                            Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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
#include "Elements/OpalMultipoleT.h"

#include "AbsBeamline/MultipoleT.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/Expressions.h"
#include "Attributes/Attributes.h"
#include "Expressions/SValue.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Units.h"

#include <iostream>
#include <sstream>
#include <vector>


OpalMultipoleT::OpalMultipoleT():
    OpalElement(SIZE, "MULTIPOLET",
                "The \"MULTIPOLET\" element defines a combined function multipole.") {
    itsAttr[TP] = Attributes::makeRealArray
        ("TP",
        "Array of multipolar field strengths b_k. The field generated in the "
        "flat top is B = b_k x^k [T m^(-k)]");

    itsAttr[LFRINGE] = Attributes::makeReal
        ("LFRINGE", "The length of the left end field [m]");

    itsAttr[RFRINGE] = Attributes::makeReal
        ("RFRINGE", "The length of the right end field [m]");

    itsAttr[HAPERT] = Attributes::makeReal
        ("HAPERT", "The aperture width [m]");

    itsAttr[VAPERT] = Attributes::makeReal
        ("VAPERT", "The aperture height [m]");

    itsAttr[ANGLE] = Attributes::makeReal
        ("ANGLE", "The azimuthal angle of the magnet in ring [rad]");

    itsAttr[EANGLE] = Attributes::makeReal
        ("EANGLE", "The entrance angle [rad]");

    itsAttr[MAXFORDER] = Attributes::makeReal
        ("MAXFORDER", "Number of terms used in each field component");

    itsAttr[MAXXORDER] = Attributes::makeReal
        ("MAXXORDER", "Number of terms used in polynomial expansions");

    itsAttr[ROTATION] = Attributes::makeReal
        ("ROTATION", "Rotation angle about its axis for skew elements [rad]");

    itsAttr[VARRADIUS] = Attributes::makeBool
        ("VARRADIUS", "Set true if radius of magnet is variable");

    itsAttr[BBLENGTH] = Attributes::makeReal
        ("BBLENGTH", "Distance between centre of magnet and entrance [m]");

    itsAttr[MAGNET_START] = Attributes::makeReal
        ("MAGNET_START", "Distance between the placement pointer and the magnet placement [m]");

    registerOwnership();
    setElement(new MultipoleT("MULTIPOLET"));
}

OpalMultipoleT::OpalMultipoleT(const std::string& name,
                               OpalMultipoleT* parent):
    OpalElement(name, parent) {
    setElement(new MultipoleT(name));
}


OpalMultipoleT::~OpalMultipoleT()
{}


OpalMultipoleT* OpalMultipoleT::clone(const std::string& name) {
    return new OpalMultipoleT(name, this);
}


void OpalMultipoleT::print(std::ostream& os) const {
    OpalElement::print(os);
}


void OpalMultipoleT::update() {
    OpalElement::update();

    MultipoleT* multT = dynamic_cast<MultipoleT*>(getElement());

    double length = Attributes::getReal(itsAttr[LENGTH]);
    double angle = Attributes::getReal(itsAttr[ANGLE]);
    multT->setElementLength(length);
    multT->setLength(length);
    multT->setBendAngle(angle);
    multT->setAperture(Attributes::getReal(itsAttr[VAPERT]),
                       Attributes::getReal(itsAttr[HAPERT]));
    multT->setFringeField(Attributes::getReal(itsAttr[LENGTH]) * 0.5,
                          Attributes::getReal(itsAttr[LFRINGE]),
                          Attributes::getReal(itsAttr[RFRINGE]));
    if (Attributes::getBool(itsAttr[VARRADIUS])) {
        multT->setVarRadius();
    }
    multT->setBoundingBoxLength(Attributes::getReal(itsAttr[BBLENGTH]));
    const std::vector<double> transProfile =  Attributes::getRealArray(itsAttr[TP]);
    int transSize = transProfile.size();

    if (transSize == 0) {
        multT->setTransMaxOrder(0);
    } else {
        multT->setTransMaxOrder(transSize - 1);
    }
    multT->setMaxOrder(Attributes::getReal(itsAttr[MAXFORDER]));
    multT->setMaxXOrder(Attributes::getReal(itsAttr[MAXXORDER]));
    multT->setRotation(Attributes::getReal(itsAttr[ROTATION]));
    multT->setEntranceAngle(Attributes::getReal(itsAttr[EANGLE]));
    multT->setMagnetStart(Attributes::getReal(itsAttr[MAGNET_START]));

    for (int comp = 0; comp < transSize; comp++) {
        double units = Units::T2kG * gsl_sf_pow_int(1.0, comp); // T m^-comp -> kG mm^-comp
        multT->setTransProfile(comp, transProfile[comp] * units);
    }
    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(multT);
    multT->initialise();

    setElement(multT);
}
