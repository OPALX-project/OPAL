//
// Class OpalVariableRFCavityFringeField
//   The class provides the user interface for
//   the VARIABLE_RF_CAVITY_FRINGE_FIELD object.
//
// Copyright (c) 2018 - 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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
#include "Elements/OpalVariableRFCavityFringeField.h"

#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "AbsBeamline/VariableRFCavityFringeField.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"


const std::string OpalVariableRFCavityFringeField::doc_string =
      std::string("The \"VARIABLE_RF_CAVITY_FRINGE_FIELD\" element defines an RF cavity ") +
      std::string("with time dependent frequency, phase and amplitude.");

OpalVariableRFCavityFringeField::OpalVariableRFCavityFringeField():
    OpalElement(SIZE, "VARIABLE_RF_CAVITY_FRINGE_FIELD", doc_string.c_str()) {
    itsAttr[PHASE_MODEL] = Attributes::makeString
        ("PHASE_MODEL",
         "The name of the phase time dependence model, which should give the phase in [rad].");

    itsAttr[AMPLITUDE_MODEL] = Attributes::makeString
        ("AMPLITUDE_MODEL",
         "The name of the amplitude time dependence model, which should give the field in [MV/m].");

    itsAttr[FREQUENCY_MODEL] = Attributes::makeString
        ("FREQUENCY_MODEL",
         "The name of the frequency time dependence model, which should give the field in [MHz].");

    itsAttr[WIDTH] = Attributes::makeReal
        ("WIDTH", "Full width of the cavity [m].");

    itsAttr[HEIGHT] = Attributes::makeReal
        ("HEIGHT", "Full height of the cavity [m].");

    itsAttr[CENTRE_LENGTH] = Attributes::makeReal
        ("CENTRE_LENGTH", "Length of the cavity field flat top [m].");

    itsAttr[END_LENGTH] = Attributes::makeReal
        ("END_LENGTH", "Length of the cavity fringe fields [m].");

    itsAttr[CAVITY_CENTRE] = Attributes::makeReal
        ("CAVITY_CENTRE", "Offset of the cavity centre from the beginning of the cavity [m].");

    itsAttr[MAX_ORDER] = Attributes::makeReal
        ("MAX_ORDER", "Maximum power of y that will be evaluated in field calculations.");

    registerOwnership();

    setElement(new VariableRFCavityFringeField("VARIABLE_RF_CAVITY_FRINGE_FIELD"));
}

OpalVariableRFCavityFringeField::OpalVariableRFCavityFringeField(
                            const std::string& name,
                            OpalVariableRFCavityFringeField* parent):
    OpalElement(name, parent) {
    VariableRFCavityFringeField* cavity = dynamic_cast<VariableRFCavityFringeField*>(parent->getElement());
    setElement(new VariableRFCavityFringeField(*cavity));
}

OpalVariableRFCavityFringeField::~OpalVariableRFCavityFringeField() {
}

OpalVariableRFCavityFringeField* OpalVariableRFCavityFringeField::clone(const std::string& name) {
    return new OpalVariableRFCavityFringeField(name, this);
}

OpalVariableRFCavityFringeField* OpalVariableRFCavityFringeField::clone() {
    return new OpalVariableRFCavityFringeField(this->getOpalName(), this);
}

void OpalVariableRFCavityFringeField::update() {
    OpalElement::update();

    VariableRFCavityFringeField* cavity = dynamic_cast<VariableRFCavityFringeField*>(getElement());

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
    cavity->setCavityCentre(cavity_centre);

    // x0 is double length of flat top so divide 2
    double centreLength = Attributes::getReal(itsAttr[CENTRE_LENGTH]);
    double endLength = Attributes::getReal(itsAttr[END_LENGTH]);
    endfieldmodel::Tanh* tanh = new endfieldmodel::Tanh(centreLength / 2.,
                                                        endLength,
                                                        maxOrder + 1);
    std::shared_ptr<endfieldmodel::EndFieldModel> end(tanh);
    cavity->setEndField(end);

    setElement(cavity);
}


size_t OpalVariableRFCavityFringeField::convertToUnsigned(double value,
                                                          const std::string& name) {
    value += unsignedTolerance; // prevent rounding error
    if (std::abs(std::floor(value) - value) > 2*unsignedTolerance) {
        throw OpalException("OpalVariableRFCavityFringeField::convertToUnsigned",
                            "Value for " + name +
                             " should be an unsigned int but a real value was found");
    }
    if (std::floor(value) < -0.5) {
        throw OpalException("OpalVariableRFCavityFringeField::convertToUnsigned",
                            "Value for " + name + " should be 0 or more");
    }
    size_t ret(std::floor(value));
    return ret;
}
