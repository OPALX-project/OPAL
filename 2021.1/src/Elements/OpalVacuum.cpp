//
// Class OpalVacuum
//   The class of OPAL Vacuum.
//
// Copyright (c) 2018-2019, Pedro Calvo, CIEMAT, Spain
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Optimizing the radioisotope production of the novel AMIT
// superconducting weak focusing cyclotron"
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
#include "Elements/OpalVacuum.h"

#include "Attributes/Attributes.h"
#include "BeamlineCore/VacuumRep.h"
#include "Physics/Physics.h"
#include "Structure/ParticleMatterInteraction.h"


OpalVacuum::OpalVacuum():
    OpalElement(SIZE, "VACUUM",
                "The \"VACUUM\" element defines the vacuum conditions for beam stripping interactions"),
    parmatint_m(NULL) {
    itsAttr[PRESSURE]        = Attributes::makeReal
        ("PRESSURE", " Pressure in the accelerator, [mbar]");
    itsAttr[TEMPERATURE]  = Attributes::makeReal
        ("TEMPERATURE", " Temperature of the accelerator, [K]");
    itsAttr[PMAPFN]   = Attributes::makeString
        ("PMAPFN", "Filename for the Pressure fieldmap");
    itsAttr[PSCALE]   = Attributes::makeReal
        ("PSCALE", "Scale factor for the P-field", 1.0);
    itsAttr[GAS]          = Attributes::makePredefinedString
        ("GAS", "The composition of residual gas", {"AIR", "H2"});
    itsAttr[STOP]         = Attributes::makeBool
        ("STOP", "Option Whether stop tracking after beam stripping. Default: true", true);
    
    registerOwnership();
    
    setElement(new VacuumRep("VACUUM"));
}


OpalVacuum::OpalVacuum(const std::string& name, OpalVacuum* parent):
    OpalElement(name, parent),
    parmatint_m(NULL) {
    setElement(new VacuumRep(name));
}


OpalVacuum::~OpalVacuum() {
    delete parmatint_m;
}


OpalVacuum* OpalVacuum::clone(const std::string& name) {
    return new OpalVacuum(name, this);
}


void OpalVacuum::update() {
    OpalElement::update();

    VacuumRep* vac =
        dynamic_cast<VacuumRep*>(getElement());

    double pressure     = Attributes::getReal(itsAttr[PRESSURE]);
    double temperature  = Attributes::getReal(itsAttr[TEMPERATURE]);
    std::string pmap    = Attributes::getString(itsAttr[PMAPFN]);
    double pscale       = Attributes::getReal(itsAttr[PSCALE]);
    bool   stop         = Attributes::getBool(itsAttr[STOP]);
    std::string gas     = Attributes::getString(itsAttr[GAS]);

    vac->setPressure(pressure);
    vac->setTemperature(temperature);
    vac->setPressureMapFN(pmap);
    vac->setPScale(pscale);
    vac->setStop(stop);
    vac->setResidualGas(gas);

    if(itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == NULL) {
        const std::string matterDescriptor = Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION]);
        ParticleMatterInteraction* orig = ParticleMatterInteraction::find(matterDescriptor);
        parmatint_m = orig->clone(matterDescriptor);
        parmatint_m->initParticleMatterInteractionHandler(*vac);
        vac->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(vac);
}
