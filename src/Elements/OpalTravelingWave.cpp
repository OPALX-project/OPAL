//
// Class OpalTravelingWave
//   The TRAVELINGWAVE element.
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
#include "Elements/OpalTravelingWave.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/TravelingWaveRep.h"
#include "Structure/OpalWake.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"

OpalTravelingWave::OpalTravelingWave():
    OpalElement(SIZE, "TRAVELINGWAVE",
                "The \"TRAVELINGWAVE\" element defines a traveling wave structure."),
    owk_m(nullptr) {
    itsAttr[VOLT] = Attributes::makeReal
                    ("VOLT", "RF voltage in MV/m");
    itsAttr[DVOLT] = Attributes::makeReal
                     ("DVOLT", "RF voltage error in MV/m");
    itsAttr[FREQ] = Attributes::makeReal
                    ("FREQ", "RF frequency in MHz");
    itsAttr[LAG] = Attributes::makeReal
                   ("LAG", "Phase lag in rad");
    itsAttr[DLAG] = Attributes::makeReal
                    ("DLAG", "Phase lag error in rad");
    itsAttr[FMAPFN] = Attributes::makeString
                      ("FMAPFN", "Filename for the fieldmap");
    itsAttr[FAST] = Attributes::makeBool
                    ("FAST", "Faster but less accurate", true);
    itsAttr[APVETO] = Attributes::makeBool
                    ("APVETO", "Do not use this cavity in the Autophase procedure", false);
    itsAttr[NUMCELLS] = Attributes::makeReal
                        ("NUMCELLS", "Number of cells in a TW structure");
    itsAttr[DESIGNENERGY] = Attributes::makeReal
                            ("DESIGNENERGY", "the mean energy of the particles at exit", -1.0);
    itsAttr[MODE] = Attributes::makeReal
                     ("MODE", "The phase shift between neighboring cells in 2*pi", 1.0/3.0);

    registerOwnership();

    setElement(new TravelingWaveRep("TRAVELINGWAVE"));
}


OpalTravelingWave::OpalTravelingWave(const std::string &name, OpalTravelingWave *parent):
    OpalElement(name, parent),
    owk_m(nullptr) {
    setElement(new TravelingWaveRep(name));
}


OpalTravelingWave::~OpalTravelingWave() {
    delete owk_m;
}


OpalTravelingWave *OpalTravelingWave::clone(const std::string &name) {
    return new OpalTravelingWave(name, this);
}


void OpalTravelingWave::update() {
    OpalElement::update();

    TravelingWaveRep *rfc =
        dynamic_cast<TravelingWaveRep *>(getElement());

    double length = Attributes::getReal(itsAttr[LENGTH]);
    double vPeak  = Attributes::getReal(itsAttr[VOLT]);
    double vPeakError  = Attributes::getReal(itsAttr[DVOLT]);
    double phase  = Attributes::getReal(itsAttr[LAG]);
    double phaseError  = Attributes::getReal(itsAttr[DLAG]);
    double freq   = Physics::two_pi * Attributes::getReal(itsAttr[FREQ]) * Units::MHz2Hz;
    std::string fmapfm = Attributes::getString(itsAttr[FMAPFN]);
    bool fast = Attributes::getBool(itsAttr[FAST]);
    bool apVeto = Attributes::getBool(itsAttr[APVETO]);

    //    std::string type = Attributes::getString(itsAttr[TYPE]);
    double kineticEnergy = Attributes::getReal(itsAttr[DESIGNENERGY]);

    rfc->setElementLength(length);
    rfc->setAmplitude(Units::MVpm2Vpm * vPeak);
    rfc->setFrequency(freq);
    rfc->setPhase(phase);

    rfc->setFieldMapFN(fmapfm);
    rfc->setFast(fast);
    rfc->setAutophaseVeto(apVeto);
    rfc->setAmplitudem(vPeak);
    rfc->setAmplitudeError(vPeakError);
    rfc->setFrequencym(freq);
    rfc->setPhasem(phase);
    rfc->setPhaseError(phaseError);
    rfc->setNumCells((int)Attributes::getReal(itsAttr[NUMCELLS]));
    rfc->setMode(Attributes::getReal(itsAttr[MODE]));
    rfc->setDesignEnergy(kineticEnergy);

    if(itsAttr[WAKEF] && owk_m == nullptr) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + std::string("_wake"));
        owk_m->initWakefunction(*rfc);
        rfc->setWake(owk_m->wf_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(rfc);
}
