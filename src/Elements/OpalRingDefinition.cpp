/*
 *  Copyright (c) 2012, Chris Rogers
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

#include "Elements/OpalRingDefinition.h"

#include <limits>

#include "AbsBeamline/Ring.h"
#include "Attributes/Attributes.h"
#include "Physics/Units.h"

OpalRingDefinition::OpalRingDefinition() :
    OpalElement(SIZE, "RINGDEFINITION",
                "The \"RINGDEFINITION\" element defines basic ring parameters.") {

    itsAttr[HARMONIC_NUMBER] = Attributes::makeReal("HARMONIC_NUMBER",
                                                    "The assumed harmonic number of the ring (i.e. number of bunches in the ring on a given turn).");
    itsAttr[LAT_RINIT] = Attributes::makeReal("LAT_RINIT",
                                                  "The initial radius of the first element to be placed in the ring [m].");
    itsAttr[LAT_PHIINIT] = Attributes::makeReal("LAT_PHIINIT", "The initial angle around the ring of the first element to be placed. [deg]");
    itsAttr[LAT_THETAINIT] = Attributes::makeReal("LAT_THETAINIT", "The angle relative to the tangent of the ring for the first element to be placed [deg].");
    itsAttr[BEAM_PHIINIT] = Attributes::makeReal("BEAM_PHIINIT",
                                                 "The initial angle around the ring of the beam [deg].");
    itsAttr[BEAM_PRINIT] = Attributes::makeReal("BEAM_PRINIT",
                                                "An initial pr momentum offset of the beam.");
    itsAttr[BEAM_RINIT] = Attributes::makeReal("BEAM_RINIT",
                                               "The initial radius of the beam [m].");
    itsAttr[SYMMETRY] = Attributes::makeReal("SYMMETRY",
                                             "The rotational symmetry of the lattice.");
    itsAttr[SCALE] = Attributes::makeReal("SCALE", "Scale the fields by a multiplicative factor",1.0);
    // should be in RF cavity definition; this comes from cyclotron definition,
    // but not right
    itsAttr[RFFREQ] = Attributes::makeReal("RFFREQ",
                                           "The nominal RF frequency of the ring [MHz].");
    // I see also makeBool, but dont know how it works; no registerBoolAttribute
    itsAttr[IS_CLOSED] = Attributes::makeBool("IS_CLOSED",
                                                "Set to 'false' to disable checking for closure of the ring");
    itsAttr[MIN_R] = Attributes::makeReal("MIN_R",
                                           "Minimum allowed radius during tracking [m]. If not defined, any radius is allowed. If MIN_R is defined, MAX_R must also be defined.");
    itsAttr[MAX_R] = Attributes::makeReal("MAX_R",
                                           "Maximum allowed radius during tracking [m]. If not defined, any radius is allowed. If MAX_R is defined, MIN_R must also be defined.");

    registerOwnership();

    setElement(new Ring("RING"));
}

OpalRingDefinition* OpalRingDefinition::clone(const std::string &name) {
    return new OpalRingDefinition(name, this);
}

void OpalRingDefinition::print(std::ostream& out) const {
    OpalElement::print(out);
}

OpalRingDefinition::OpalRingDefinition(const std::string &name, OpalRingDefinition *parent):
    OpalElement(name, parent) {
    setElement(new Ring(name));
}

OpalRingDefinition::~OpalRingDefinition() {}

void OpalRingDefinition::update() {
    Ring *ring = dynamic_cast<Ring*>(getElement());
    ring->setBeamPhiInit(Attributes::getReal(itsAttr[BEAM_PHIINIT]));
    ring->setBeamPRInit(Attributes::getReal(itsAttr[BEAM_PRINIT]));
    ring->setBeamRInit(Attributes::getReal(itsAttr[BEAM_RINIT])*Units::m2mm);
    ring->setLatticeRInit(Attributes::getReal(itsAttr[LAT_RINIT])*Units::m2mm);

    ring->setLatticePhiInit(Attributes::getReal(itsAttr[LAT_PHIINIT])*Units::deg2rad);
    ring->setLatticeThetaInit(Attributes::getReal(itsAttr[LAT_THETAINIT])*Units::deg2rad);
    ring->setSymmetry(Attributes::getReal(itsAttr[SYMMETRY]));
    ring->setScale(Attributes::getReal(itsAttr[SCALE]));

    ring->setHarmonicNumber(Attributes::getReal(itsAttr[HARMONIC_NUMBER]));
    ring->setRFFreq(Attributes::getReal(itsAttr[RFFREQ]));
    ring->setIsClosed(Attributes::getBool(itsAttr[IS_CLOSED]));
    double minR = -1;
    double maxR = -1;

    if (itsAttr[MIN_R]) {
        minR = Attributes::getReal(itsAttr[MIN_R]);
        if (!itsAttr[MAX_R]) {
            throw (""); // EXCEPTION
        }
    }
    if (itsAttr[MAX_R]) {
        maxR = Attributes::getReal(itsAttr[MAX_R]);
        if (!itsAttr[MIN_R]) {
            throw (""); // EXCEPTION
        }
        ring->setRingAperture(minR, maxR);
    }

    setElement(ring);
}