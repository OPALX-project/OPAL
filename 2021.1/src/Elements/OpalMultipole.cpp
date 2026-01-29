//
// Class OpalMultipole
//   The MULTIPOLE element.
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
#include "Elements/OpalMultipole.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MultipoleRep.h"
#include "Expressions/SValue.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include <iostream>
#include <sstream>
#include <vector>

OpalMultipole::OpalMultipole():
    OpalElement(SIZE, "MULTIPOLE",
                "The \"MULTIPOLE\" element defines a thick multipole.\n"
                "* If the length is non-zero, the strengths are per unit "
                "length.\n* If the length is zero, the strengths are the "
                "values integrated over the length.\n"
                "* With zero length no synchrotron radiation can be calculated.") {
    itsAttr[KN] = Attributes::makeRealArray
                  ("KN", "Normalised multipole strengths (normal) in m^(-k)");
    itsAttr[DKN] = Attributes::makeRealArray
                  ("DKN", "Normalised multipole strengths errors(normal) in m^(-k)");
    itsAttr[KS] = Attributes::makeRealArray
                  ("KS", "Normalised multipole strengths (skew) in m^(-k)");
    itsAttr[DKS] = Attributes::makeRealArray
                  ("DKS", "Normalised multipole strength errors (skew) in m^(-k)");

    registerOwnership();

    setElement(new MultipoleRep("MULTIPOLE"));
}


OpalMultipole::OpalMultipole(const std::string &name, OpalMultipole *parent):
    OpalElement(name, parent) {
    setElement(new MultipoleRep(name));
}


OpalMultipole::~OpalMultipole()
{}


OpalMultipole *OpalMultipole::clone(const std::string &name) {
    return new OpalMultipole(name, this);
}


void OpalMultipole::print(std::ostream &os) const {
    OpalElement::print(os);
}


void OpalMultipole::update() {
    OpalElement::update();

    // Magnet length.
    MultipoleRep *mult =
        dynamic_cast<MultipoleRep *>(getElement());
    double length = getLength();
    mult->setElementLength(length);

    // Field components.
    BMultipoleField field;

    const std::vector<double> norm = Attributes::getRealArray(itsAttr[KN]);
    std::vector<double> normErrors = Attributes::getRealArray(itsAttr[DKN]);
    const std::vector<double> skew = Attributes::getRealArray(itsAttr[KS]);
    std::vector<double> skewErrors = Attributes::getRealArray(itsAttr[DKS]);
    int normSize = norm.size();
    int skewSize = skew.size();
    normErrors.resize(normSize, 0.0);
    skewErrors.resize(skewSize, 0.0);
    double factor = OpalData::getInstance()->getP0() / Physics::c;
    int top = (normSize > skewSize) ? normSize : skewSize;

    for(int comp = 1; comp <= top; comp++) {
        factor /= double(comp);
        if(comp <= normSize) {
            field.setNormalComponent(comp, norm[comp-1] * factor);
            mult->setNormalComponent(comp, norm[comp-1], normErrors[comp-1]);
        }
        if(comp <= skewSize) {
            field.setSkewComponent(comp, skew[comp-1] * factor);
            mult->setSkewComponent(comp, skew[comp-1], skewErrors[comp-1]);
        }
    }

    mult->setField(field);

   // Transmit "unknown" attributes.
    OpalElement::updateUnknown(mult);
}