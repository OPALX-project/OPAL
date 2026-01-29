/*
 *  Copyright (c) 2017, Titus Dascalu
 *  Copyright (c) 2018, Martin Duy Tat
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */


#include "Elements/OpalMultipoleTStraight.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Expressions/SValue.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include <iostream>
#include <sstream>
#include <vector>


// Class OpalMultipoleTStraight
// ------------------------------------------------------------------------

OpalMultipoleTStraight::OpalMultipoleTStraight():
    OpalElement(SIZE, "MULTIPOLETSTRAIGHT",
    "The \"MULTIPOLETSTRAIGHT\" element defines a straight, combined function multipole magnet.") {
    itsAttr[TP] = Attributes::makeRealArray
                  ("TP", "Transverse Profile derivatives in m^(-k)");
    itsAttr[LFRINGE] = Attributes::makeReal
                  ("LFRINGE", "The length of the left end field in m");
    itsAttr[RFRINGE] = Attributes::makeReal
                  ("RFRINGE", "The length of the right end field in m");
    itsAttr[HAPERT] = Attributes::makeReal
                  ("HAPERT", "The aperture width in m");
    itsAttr[VAPERT] = Attributes::makeReal
                  ("VAPERT", "The aperture height in m");
    itsAttr[EANGLE] = Attributes::makeReal
                  ("EANGLE", "The entrance angle (rad)");
    itsAttr[MAXFORDER] = Attributes::makeReal
                  ("MAXFORDER",
                   "Number of terms used in each field component");
    itsAttr[ROTATION] = Attributes::makeReal
                  ("ROTATION",
                   "Rotation angle about its axis for skew elements (rad)");
    itsAttr[BBLENGTH] = Attributes::makeReal
                  ("BBLENGTH",
                   "Distance between centre of magnet and entrance in m");

    registerOwnership();

    setElement(new MultipoleTStraight("MULTIPOLETSTRAIGHT"));
}


OpalMultipoleTStraight::OpalMultipoleTStraight(const std::string &name,
                                               OpalMultipoleTStraight *parent):
    OpalElement(name, parent) {
    setElement(new MultipoleTStraight(name));
}


OpalMultipoleTStraight::~OpalMultipoleTStraight()
{}


OpalMultipoleTStraight *OpalMultipoleTStraight::clone(const std::string &name) {
    return new OpalMultipoleTStraight(name, this);
}


void OpalMultipoleTStraight::print(std::ostream &os) const {
    OpalElement::print(os);
}


void OpalMultipoleTStraight::update() {
    OpalElement::update();

    // Magnet length.
    MultipoleTStraight *multT =
    dynamic_cast<MultipoleTStraight*>(getElement());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double boundingBoxLength = Attributes::getReal(itsAttr[BBLENGTH]);
    multT->setElementLength(length);
    multT->setLength(length);
    multT->setAperture(Attributes::getReal(itsAttr[VAPERT]),
                       Attributes::getReal(itsAttr[HAPERT]));

    multT->setFringeField(Attributes::getReal(itsAttr[LENGTH])/2,
                          Attributes::getReal(itsAttr[LFRINGE]),
                          Attributes::getReal(itsAttr[RFRINGE]));
    multT->setBoundingBoxLength(Attributes::getReal(itsAttr[BBLENGTH]));
    const std::vector<double> transProfile =
                              Attributes::getRealArray(itsAttr[TP]);
    int transSize = transProfile.size();

    multT->setTransMaxOrder(transSize - 1);
    multT->setMaxOrder(Attributes::getReal(itsAttr[MAXFORDER]));
    multT->setRotation(Attributes::getReal(itsAttr[ROTATION]));
    multT->setEntranceAngle(Attributes::getReal(itsAttr[EANGLE]));

    StraightGeometry &geometry = multT->getGeometry();

    geometry = StraightGeometry(2 * boundingBoxLength);

    for(int comp = 0; comp < transSize; comp++) {
        multT->setTransProfile(comp, transProfile[comp]);
    }
    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(multT);

    setElement(multT);
}