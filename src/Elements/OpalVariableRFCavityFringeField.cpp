/*
 *  Copyright (c) 2018, Chris Rogers
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

#include "Elements/OpalVariableRFCavityFringeField.h"

#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "AbsBeamline/VariableRFCavityFringeField.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "Attributes/Attributes.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Utilities/OpalException.h"
#include "Utilities/OpalException.h"


const std::string OpalVariableRFCavityFringeField::doc_string =
      std::string("The \"VARIABLE_RF_CAVITY_FRINGE_FIELD\" element defines an RF cavity ")+
      std::string("with time dependent frequency, phase and amplitude.");

OpalVariableRFCavityFringeField::OpalVariableRFCavityFringeField():
    OpalElement(SIZE, "VARIABLE_RF_CAVITY_FRINGE_FIELD", doc_string.c_str()) {
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
    itsAttr[CENTRE_LENGTH] = Attributes::makeReal("CENTRE_LENGTH",
                "Length of the cavity field flat top [m].");
    itsAttr[END_LENGTH] = Attributes::makeReal("END_LENGTH",
                "Length of the cavity fringe fields [m].");
    itsAttr[CAVITY_CENTRE] = Attributes::makeReal("CAVITY_CENTRE",
                "Offset of the cavity centre from the beginning of the cavity [m].");
    itsAttr[MAX_ORDER] = Attributes::makeReal("MAX_ORDER",
                "Maximum power of y that will be evaluated in field calculations.");

    registerOwnership();

    setElement(new VariableRFCavityFringeField("VARIABLE_RF_CAVITY_FRINGE_FIELD"));
}

OpalVariableRFCavityFringeField::OpalVariableRFCavityFringeField(
                            const std::string &name,
                            OpalVariableRFCavityFringeField *parent
                        ) : OpalElement(name, parent) {
    VariableRFCavityFringeField *cavity = dynamic_cast
         <VariableRFCavityFringeField*>(parent->getElement());
    setElement(new VariableRFCavityFringeField(*cavity));
}

OpalVariableRFCavityFringeField::~OpalVariableRFCavityFringeField() {
}

OpalVariableRFCavityFringeField *OpalVariableRFCavityFringeField::clone(
                                                      const std::string &name) {
    return new OpalVariableRFCavityFringeField(name, this);
}

OpalVariableRFCavityFringeField *OpalVariableRFCavityFringeField::clone() {
    return new OpalVariableRFCavityFringeField(this->getOpalName(), this);
}

void OpalVariableRFCavityFringeField::update() {
    OpalElement::update();

    VariableRFCavityFringeField *cavity = dynamic_cast
                <VariableRFCavityFringeField*>(getElement());
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
    double maxOrderReal = Attributes::getReal(itsAttr[MAX_ORDER]);
    size_t maxOrder = convertToUnsigned(maxOrderReal, "MAX_ORDER");
    cavity->setMaxOrder(maxOrder);
    double cavity_centre = Attributes::getReal(itsAttr[CAVITY_CENTRE]);
    cavity->setCavityCentre(cavity_centre); // mm
    // convert to mm; x0 is double length of flat top so divide 2
    double centreLength = Attributes::getReal(itsAttr[CENTRE_LENGTH]) * Units::m2mm;
    double endLength = Attributes::getReal(itsAttr[END_LENGTH]) * Units::m2mm;
    endfieldmodel::Tanh* tanh = new endfieldmodel::Tanh(centreLength/2.,
                                                        endLength,
                                                        maxOrder+1);
    std::shared_ptr<endfieldmodel::EndFieldModel> end(tanh);
    cavity->setEndField(end);

    setElement(cavity);
}


size_t OpalVariableRFCavityFringeField::convertToUnsigned(double value,
                                                          std::string name) {
    value += unsignedTolerance; // prevent rounding error
    if (std::abs(std::floor(value) - value) > 2*unsignedTolerance) {
        throw OpalException("OpalVariableRFCavityFringeField::convertToUnsigned",
                    "Value for "+name+
                    " should be an unsigned int but a real value was found");
    }
    if (std::floor(value) < -0.5) {
        throw OpalException("OpalVariableRFCavityFringeField::convertToUnsigned",
                            "Value for "+name+" should be 0 or more");
    }
    size_t ret(std::floor(value));
    return ret;
}