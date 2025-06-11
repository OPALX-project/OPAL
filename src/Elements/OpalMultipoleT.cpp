//
// Class OpalMultipoleT
//   The Opal MultipoleT element.
//
// Copyright (c) 2017 - 2023, Titus Dascalu
//                            Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
// Copyright (c) 2025, Jon Thompson
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
#include <iostream>
#include <vector>
#include "AbsBeamline/MultipoleT.h"
#include "Attributes/Attributes.h"

OpalMultipoleT::OpalMultipoleT()
    : OpalElement(
          SIZE, "MULTIPOLET", "The \"MULTIPOLET\" element defines a combined function multipole.") {
    // Attributes for a straight multipole
    itsAttr[TP] = Attributes::makeRealArray(
        "TP",
        "Array of multipolar field strengths b_k. The field generated in the "
        "flat top is B = b_k x^k [T m^(-k)]");
    itsAttr[LFRINGE] = Attributes::makeReal("LFRINGE", "The length of the left end field [m]");
    itsAttr[RFRINGE] = Attributes::makeReal("RFRINGE", "The length of the right end field [m]");
    itsAttr[HAPERT]  = Attributes::makeReal("HAPERT", "The aperture width [m]");
    itsAttr[VAPERT]  = Attributes::makeReal("VAPERT", "The aperture height [m]");
    itsAttr[MAXFORDER] = Attributes::makeReal(
        "MAXFORDER", "Number of terms used in each fringe component", DefaultMAXFORDER);
    itsAttr[ROTATION] = Attributes::makeReal(
        "ROTATION", "Rotation angle about its axis for skew elements [rad]");
    itsAttr[EANGLE] = Attributes::makeReal("EANGLE", "The entrance angle [rad]");
    itsAttr[BBLENGTH] = Attributes::makeReal(
        "BBLENGTH", "Distance between centre of magnet and entrance [m]");

    // Further attributes for a constant radius curved multipole
    itsAttr[ANGLE] = Attributes::makeReal(
        "ANGLE", "The azimuthal angle of the magnet in ring [rad]", 0.0);
    itsAttr[MAXXORDER] = Attributes::makeReal(
        "MAXXORDER", "Number of terms used in polynomial expansions", DefaultMAXXORDER);

    // Further attributes for a variable radius multipole
    itsAttr[VARRADIUS] = Attributes::makeBool(
        "VARRADIUS", "Set true if radius of magnet is variable", false);
    itsAttr[ENTRYOFFSET] = Attributes::makeReal(
        "ENTRYOFFSET", "Longitudinal offset from standard entrance point [m]", 0.0);

    registerOwnership();
    setElement(new MultipoleT("MULTIPOLET"));
}

OpalMultipoleT::OpalMultipoleT(const std::string& name, OpalMultipoleT* parent)
    : OpalElement(name, parent) {
    setElement(new MultipoleT(name));
}

OpalMultipoleT* OpalMultipoleT::clone(const std::string& name) {
    return new OpalMultipoleT(name, this);
}

void OpalMultipoleT::print(std::ostream& os) const {
    OpalElement::print(os);
}

void OpalMultipoleT::update() {
    // Base class first
    OpalElement::update();
    // Make some sanity checks
    auto maxFOrder = Attributes::getReal(itsAttr[MAXFORDER]);
    if(maxFOrder < MinimumMAXFORDER) {
        throw OpalException("OpalMultipoleT::Update",
                            "Attribute MAXFORDER must be >= 1.0");
    }
    if(maxFOrder > MaximumMAXFORDER) {
        WARNMSG("OpalMultipoleT::Update, a value of "
                << maxFOrder << " for MAXFORDER may lead to excessive run time");
    }
    auto rotation = Attributes::getReal(itsAttr[ROTATION]);
    double bendAngle = Attributes::getReal(itsAttr[ANGLE]);
    if(bendAngle != 0.0 && rotation != 0.0) {
        throw OpalException("OpalMultipoleT::Update",
                            "Non-zero ROTATION (a skew multipole) is only supported for straight magnets");
    }
    bool varRadius = Attributes::getBool(itsAttr[VARRADIUS]);
    if(varRadius && bendAngle != 0.0) {
        WARNMSG("OpalMultipoleT::Update, the variable radius multipole magnet implementation is very slow");
    }
    double entryOffset = Attributes::getReal(itsAttr[ENTRYOFFSET]);
    if((!varRadius || bendAngle == 0.0) && entryOffset != 0.0) {
        throw OpalException("OpalMultipoleT::Update",
                            "The ENTRYOFFSET is only supported for variable radius curved magnets");
    }
    // Convert pole strengths from Tesla to internal units which are kGauss
    auto tp = Attributes::getRealArray(itsAttr[TP]);
    for(auto& i : tp) {
        i *= Units::T2kG;
    }
    // Set the attributes
    auto length = Attributes::getReal(itsAttr[LENGTH]);
    auto* multT = dynamic_cast<MultipoleT*>(getElement());
    multT->setElementLength(length);
    multT->setBendAngle(bendAngle, varRadius);
    multT->setAperture(Attributes::getReal(itsAttr[VAPERT]), Attributes::getReal(itsAttr[HAPERT]));
    multT->setFringeField(
        length * 0.5, Attributes::getReal(itsAttr[LFRINGE]), Attributes::getReal(itsAttr[RFRINGE]));
    multT->setBoundingBoxLength(Attributes::getReal(itsAttr[BBLENGTH]));
    multT->setTransProfile(tp);
    multT->setMaxOrder(static_cast<size_t>(maxFOrder),
        static_cast<size_t>(Attributes::getReal(itsAttr[MAXXORDER])));
    multT->setRotation(rotation);
    multT->setEntranceAngle(Attributes::getReal(itsAttr[EANGLE]));
    multT->setEntryOffset(entryOffset);
    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(multT);
}
