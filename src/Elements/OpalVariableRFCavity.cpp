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

#include "Elements/OpalVariableRFCavity.h"

#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "Attributes/Attributes.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "AbsBeamline/VariableRFCavity.h"

extern Inform *gmsg;

const std::string OpalVariableRFCavity::doc_string =
      std::string("The \"VARIABLE_RF_CAVITY\" element defines an RF cavity ")+
      std::string("with time dependent frequency, phase and amplitude.");

OpalVariableRFCavity::OpalVariableRFCavity():
    OpalElement(SIZE, "VARIABLE_RF_CAVITY", doc_string.c_str()) {
    itsAttr[PHASE_MODEL] = Attributes::makeString("PHASE_MODEL",
                "The name of the phase time dependence model, which should give the phase in [rad].");
    itsAttr[AMPLITUDE_MODEL] = Attributes::makeString("AMPLITUDE_MODEL",
                "The name of the amplitude time dependence model, which should give the field in [MV/m]");
    itsAttr[FREQUENCY_MODEL] = Attributes::makeString("FREQUENCY_MODEL",
                "The name of the frequency time dependence model, which should give the field in [MHz].");
    itsAttr[WIDTH] = Attributes::makeReal("WIDTH",
                "Full width of the cavity [m].");
    itsAttr[HEIGHT] = Attributes::makeReal("HEIGHT",
                "Full height of the cavity [m].");

    registerOwnership();

    setElement(new VariableRFCavity("VARIABLE_RF_CAVITY"));
}

OpalVariableRFCavity::OpalVariableRFCavity(const std::string &name,
                                           OpalVariableRFCavity *parent) :
          OpalElement(name, parent) {
    VariableRFCavity *cavity = dynamic_cast<VariableRFCavity*>(
                                        parent->getElement());
    setElement(new VariableRFCavity(*cavity));
}

OpalVariableRFCavity::~OpalVariableRFCavity() {
}

OpalVariableRFCavity *OpalVariableRFCavity::clone(const std::string &name) {
    return new OpalVariableRFCavity(name, this);
}

OpalVariableRFCavity *OpalVariableRFCavity::clone() {
    return new OpalVariableRFCavity(this->getOpalName(), this);
}

void OpalVariableRFCavity::update() {
    OpalElement::update();

    VariableRFCavity *cavity = dynamic_cast<VariableRFCavity*>(
                                                getElement());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    cavity->setLength(length);
    std::string phaseName = Attributes::getString(itsAttr[PHASE_MODEL]);
    cavity->setPhaseName(phaseName);
    std::string ampName = Attributes::getString(itsAttr[AMPLITUDE_MODEL]);
    cavity->setAmplitudeName(ampName);
    std::string freqName = Attributes::getString(itsAttr[FREQUENCY_MODEL]);
    cavity->setFrequencyName(freqName);
    double width = Attributes::getReal(itsAttr[WIDTH]);
    cavity->setWidth(width);
    double height = Attributes::getReal(itsAttr[HEIGHT]);
    cavity->setHeight(height);
    setElement(cavity);
}