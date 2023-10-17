//
// Class OpalVariableRFCavity
//   The class provides the user interface for the VARIABLE_RF_CAVITY object.
//
// Copyright (c) 2014 - 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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
#include "Elements/OpalVariableRFCavity.h"

#include "AbsBeamline/VariableRFCavity.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "Attributes/Attributes.h"

const std::string OpalVariableRFCavity::doc_string =
      std::string("The \"VARIABLE_RF_CAVITY\" element defines an RF cavity ") +
      std::string("with time dependent frequency, phase and amplitude.");

OpalVariableRFCavity::OpalVariableRFCavity():
    OpalElement(SIZE, "VARIABLE_RF_CAVITY", doc_string.c_str()) {

    itsAttr[PHASE_MODEL] = Attributes::makeString
        ("PHASE_MODEL",
         "The name of the phase time dependence model, which should give the phase in [rad].");

    itsAttr[AMPLITUDE_MODEL] = Attributes::makeString
        ("AMPLITUDE_MODEL",
         "The name of the amplitude time dependence model, which should give the field in [MV/m]");

    itsAttr[FREQUENCY_MODEL] = Attributes::makeString
        ("FREQUENCY_MODEL",
         "The name of the frequency time dependence model, which should give the field in [MHz].");

    itsAttr[WIDTH] = Attributes::makeReal
        ("WIDTH",
         "Full width of the cavity [m].");

    itsAttr[HEIGHT] = Attributes::makeReal
        ("HEIGHT",
         "Full height of the cavity [m].");

    registerOwnership();

    setElement(new VariableRFCavity("VARIABLE_RF_CAVITY"));
}

OpalVariableRFCavity::OpalVariableRFCavity(const std::string& name,
                                           OpalVariableRFCavity* parent):
          OpalElement(name, parent) {
    VariableRFCavity* cavity = dynamic_cast<VariableRFCavity*>(parent->getElement());
    setElement(new VariableRFCavity(*cavity));
}

OpalVariableRFCavity::~OpalVariableRFCavity() {
}

OpalVariableRFCavity* OpalVariableRFCavity::clone(const std::string& name) {
    return new OpalVariableRFCavity(name, this);
}

OpalVariableRFCavity* OpalVariableRFCavity::clone() {
    return new OpalVariableRFCavity(this->getOpalName(), this);
}

void OpalVariableRFCavity::update() {
    OpalElement::update();

    VariableRFCavity*cavity = dynamic_cast<VariableRFCavity*>(getElement());

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
