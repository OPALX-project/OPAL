//
// Class OpalUndulator
//   Defines the Undulator/Wiggler element and its attributes.
//
// Copyright (c) 2020, Arnau Albà, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// Implemented as part of the MSc thesis
// "Start-to-End Modelling of the AWA Micro-Bunched Electron Cooling POP-Experiment"
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
#include "Elements/OpalUndulator.h"

#include "Attributes/Attributes.h"
#include "BeamlineCore/UndulatorRep.h"

OpalUndulator::OpalUndulator()
    : OpalElement(SIZE, "UNDULATOR", "The \"UNDULATOR\" element defines an undulator or wiggler.") {
    itsAttr[K] = Attributes::makeReal("K", "Undulator strength parameter [unitless].", 1);

    itsAttr[LAMBDA] = Attributes::makeReal("LAMBDA", "Undulator period [m].", 0.0);

    itsAttr[NUMPERIODS] = Attributes::makeReal("NUMPERIODS", "Number of undulator periods.", 0.0);

    itsAttr[ANGLE] = Attributes::makeReal(
        "ANGLE", "Polarisation angle of the undulator magnetic field [rad].", 0.0);

    itsAttr[FNAME] = Attributes::makeString(
        "FNAME", "Jobfile specifying the output data from the undulator.", "");

    itsAttr[MESHLENGTH] =
        Attributes::makeRealArray("MESHLENGTH", "Size of computational mesh [m].");

    itsAttr[MESHRESOLUTION] =
        Attributes::makeRealArray("MESHRESOLUTION", "{dx, dy, dz} of the mesh [m].");

    itsAttr[TRUNORDER] =
        Attributes::makeReal("TRUNORDER", "Order of absorbing boundary conditions. 1st or 2nd.", 2);

    itsAttr[TOTALTIME] =
        Attributes::makeReal("TOTALTIME", "Total time of the full-wave simulation [s].", 0.0);

    itsAttr[DTBUNCH] = Attributes::makeReal(
        "DTBUNCH",
        "Time step for the particle update can be smaller than the field update step [s].", 0.0);

    registerOwnership();

    setElement(new UndulatorRep("UNDULATOR"));
}

OpalUndulator::OpalUndulator(const std::string& name, OpalUndulator* parent)
    : OpalElement(name, parent) {
    setElement(new UndulatorRep(name));
}

OpalUndulator::~OpalUndulator() {
}

OpalUndulator* OpalUndulator::clone(const std::string& name) {
    return new OpalUndulator(name, this);
}

void OpalUndulator::update() {
    OpalElement::update();

    UndulatorRep* ur = static_cast<UndulatorRep*>(getElement());
    /* The element length is given by the number and length of the undulator periods, plus the
     * length of the two fringe fields, each measuring 2 * lambda. */
    ur->setElementLength(
        Attributes::getReal(itsAttr[LAMBDA]) * (4 + Attributes::getReal(itsAttr[NUMPERIODS])));
    ur->setK(Attributes::getReal(itsAttr[K]));
    ur->setLambda(Attributes::getReal(itsAttr[LAMBDA]));
    ur->setNumPeriods(Attributes::getReal(itsAttr[NUMPERIODS]));
    ur->setAngle(Attributes::getReal(itsAttr[ANGLE]));
    ur->setFilename(Attributes::getString(itsAttr[FNAME]));
    ur->setMeshLength(Attributes::getRealArray(itsAttr[MESHLENGTH]));
    ur->setMeshResolution(Attributes::getRealArray(itsAttr[MESHRESOLUTION]));
    ur->setTruncationOrder(Attributes::getReal(itsAttr[TRUNORDER]));
    ur->setTotalTime(Attributes::getReal(itsAttr[TOTALTIME]));
    ur->setDtBunch(Attributes::getReal(itsAttr[DTBUNCH]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(ur);
}
