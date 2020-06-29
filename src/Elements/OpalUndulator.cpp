//
// Class OpalUndulator
// Defines the Undulator/Wiggler element and its attributes.
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

OpalUndulator::OpalUndulator():
    OpalElement(SIZE, "UNDULATOR",
                "The \"UNDULATOR\" element defines a undulator.") {

    itsAttr[NSLICES] = Attributes::makeReal
                          ("NSLICES",
                          "The number of slices/ steps for this element in Map Tracking.", 1);

    itsAttr[K] = Attributes::makeReal
                          ("K",
                          "The undulator parameter.", 1);

    itsAttr[LAMBDA] = Attributes::makeReal
                          ("LAMBDA",
                          "The undulator period.", 0.0);
    
    itsAttr[NUMPERIODS] = Attributes::makeReal
                          ("NUMPERIODS",
                          "Number of undulator period.", 0.0);

    itsAttr[ANGLE] = Attributes::makeReal
                          ("ANGLE",
                          "Polarisation angle of the undulator magnetic field.", 0.0);
    
    itsAttr[FNAME] = Attributes::makeString
                          ("FNAME",
                           "Jobfile specifying the output data from the undulator.", "");
    
    itsAttr[MESHLENGTH] = Attributes::makeRealArray
                          ("MESHLENGTH",
                           "Size of computational mesh.");
    
    itsAttr[MESHRESOLUTION] = Attributes::makeRealArray
                          ("MESHRESOLUTION",
                           "dx, dy, dz of the mesh.");
    
    itsAttr[TRUNORDER] = Attributes::makeReal
                          ("TRUNORDER",
                          "Order of boundary absorbing conditions. 1st or 2nd.", 2);
    
    itsAttr[TOTALTIME] = Attributes::makeReal
                          ("TOTALTIME",
                          "Total time of undulator simulation.", 0.0);
    
    itsAttr[DTBUNCH] = Attributes::makeReal
                          ("DTBUNCH",
                          "Time step for the bunch position update can be smaller than field update step.", 0.0);

    registerRealAttribute("NSLICES");
    registerRealAttribute("K");
    registerRealAttribute("LAMBDA");
    registerRealAttribute("NUMPERIODS");
    registerRealAttribute("ANGLE");
    registerStringAttribute("FNAME");
    registerRealAttribute("TRUNORDER");
    registerRealAttribute("TOTALTIME");
    registerRealAttribute("DTBUNCH");
    registerOwnership();

    setElement(new UndulatorRep("UNDULATOR"));
}


OpalUndulator::OpalUndulator(const std::string &name, OpalUndulator *parent):
    OpalElement(name, parent) {
    setElement(new UndulatorRep(name));
}


OpalUndulator::~OpalUndulator()
{}


OpalUndulator *OpalUndulator::clone(const std::string &name) {
    return new OpalUndulator(name, this);
}


void OpalUndulator::update() {
    OpalElement::update();

    UndulatorRep *ur = static_cast<UndulatorRep *>(getElement());
    ur->setElementLength(Attributes::getReal(itsAttr[LAMBDA]) * (4 + Attributes::getReal(itsAttr[NUMPERIODS])));
    ur->setNSlices(Attributes::getReal(itsAttr[NSLICES]));
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
