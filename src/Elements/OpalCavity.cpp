//
// Class OpalCavity
//   The RFCAVITY element.
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
#include "Elements/OpalCavity.h"

#include "AbstractObjects/Attribute.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/RFCavityRep.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/OpalWake.h"


OpalCavity::OpalCavity():
    OpalElement(SIZE, "RFCAVITY",
                "The \"RFCAVITY\" element defines an RF cavity."),
    owk_m(nullptr),
    obgeo_m(nullptr) {
    itsAttr[VOLT] = Attributes::makeReal
                    ("VOLT", "RF voltage in MV");
    itsAttr[DVOLT] = Attributes::makeReal
                     ("DVOLT", "RF voltage error in MV");
    itsAttr[FREQ] = Attributes::makeReal
                    ("FREQ", "RF frequency in MHz");
    itsAttr[LAG] = Attributes::makeReal
                   ("LAG", "Phase lag (rad)");
    itsAttr[DLAG] = Attributes::makeReal
                    ("DLAG", "Phase lag error (rad)");
    itsAttr[FMAPFN] = Attributes::makeString
                      ("FMAPFN", "Filename of the fieldmap");
    itsAttr[GEOMETRY] = Attributes::makeString
                        ("GEOMETRY", "BoundaryGeometry for Cavities");
    itsAttr[FAST] = Attributes::makeBool
                    ("FAST", "Faster but less accurate", true);
    itsAttr[APVETO] = Attributes::makeBool
                    ("APVETO", "Do not use this cavity in the Autophase procedure", false);
    itsAttr[RMIN] = Attributes::makeReal
                    ("RMIN", " Minimal Radius of a cyclotron cavity [mm]");
    itsAttr[RMAX] = Attributes::makeReal
                    ("RMAX", " Maximal Radius of a cyclotron cavity [mm]");
    itsAttr[ANGLE] = Attributes::makeReal
                     ("ANGLE", "Azimuth position of a cyclotron cavity [deg]");
    itsAttr[PDIS] = Attributes::makeReal
                    ("PDIS", "Shift distance of cavity gap from center of cyclotron [mm]");
    itsAttr[GAPWIDTH] = Attributes::makeReal
                        ("GAPWIDTH", "Gap width of a cyclotron cavity [mm]");
    itsAttr[PHI0] = Attributes::makeReal
                    ("PHI0", "Initial phase of cavity [deg]");
    itsAttr[DESIGNENERGY] = Attributes::makeReal
                            ("DESIGNENERGY", "the mean energy of the particles at exit", -1.0);
    // attibutes for timedependent values
    itsAttr[PHASE_MODEL] = Attributes::makeString("PHASE_MODEL",
                                                  "The name of the phase time dependence model.");
    itsAttr[AMPLITUDE_MODEL] = Attributes::makeString("AMPLITUDE_MODEL",
                                                      "The name of the amplitude time dependence model.");
    itsAttr[FREQUENCY_MODEL] = Attributes::makeString("FREQUENCY_MODEL",
                                                      "The name of the frequency time dependence model.");

    registerOwnership();

    setElement(new RFCavityRep("RFCAVITY"));
}


OpalCavity::OpalCavity(const std::string &name, OpalCavity *parent):
    OpalElement(name, parent),
    owk_m(nullptr),
    obgeo_m(nullptr) {
    setElement(new RFCavityRep(name));
}


OpalCavity::~OpalCavity() {
    delete owk_m;
}


OpalCavity *OpalCavity::clone(const std::string &name) {
    return new OpalCavity(name, this);
}


void OpalCavity::update() {
    OpalElement::update();

    RFCavityRep *rfc =
        dynamic_cast<RFCavityRep *>(getElement());

    double length = Attributes::getReal(itsAttr[LENGTH]);
    double peak = Attributes::getReal(itsAttr[VOLT]);
    double peakError = Attributes::getReal(itsAttr[DVOLT]);
    double phase = Attributes::getReal(itsAttr[LAG]);
    double phaseError = Attributes::getReal(itsAttr[DLAG]);
    double freq = Physics::two_pi * Attributes::getReal(itsAttr[FREQ]) * Units::MHz2Hz;
    std::string fmapfn = Attributes::getString(itsAttr[FMAPFN]);
    std::string type = Attributes::getString(itsAttr[TYPE]);
    bool fast = Attributes::getBool(itsAttr[FAST]);
    bool apVeto = (Attributes::getBool(itsAttr[APVETO]));

    double rmin = Attributes::getReal(itsAttr[RMIN]);
    double rmax = Attributes::getReal(itsAttr[RMAX]);
    double angle = Attributes::getReal(itsAttr[ANGLE]);
    double pdis = Attributes::getReal(itsAttr[PDIS]);
    double gapwidth = Attributes::getReal(itsAttr[GAPWIDTH]);
    double phi0 = Attributes::getReal(itsAttr[PHI0]);
    double kineticEnergy = Attributes::getReal(itsAttr[DESIGNENERGY]);

    if(itsAttr[WAKEF] && owk_m == nullptr) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + std::string("_wake"));
        owk_m->initWakefunction(*rfc);
        rfc->setWake(owk_m->wf_m);
    }

    if(itsAttr[GEOMETRY] && obgeo_m == nullptr) {
        obgeo_m = (BoundaryGeometry::find(Attributes::getString(itsAttr[GEOMETRY])))->clone(getOpalName() + std::string("_geometry"));
        if(obgeo_m) {
            rfc->setBoundaryGeometry(obgeo_m);
        }
    }

    rfc->setElementLength(length);

    rfc->setAmplitude(Units::MVpm2Vpm * peak);
    rfc->setFrequency(freq);
    rfc->setPhase(phase);

    rfc->dropFieldmaps();

    rfc->setAmplitudem(peak);
    rfc->setAmplitudeError(peakError);
    rfc->setFrequencym(freq);
    rfc->setPhasem(phase);
    rfc->setPhaseError(phaseError);
    rfc->setFieldMapFN(fmapfn);

    rfc->setFast(fast);
    rfc->setAutophaseVeto(apVeto);
    rfc->setCavityType(type);
    rfc->setRmin(rmin);
    rfc->setRmax(rmax);
    rfc->setAzimuth(angle);
    rfc->setPerpenDistance(pdis);
    rfc->setGapWidth(gapwidth);
    rfc->setPhi0(phi0);
    rfc->setDesignEnergy(kineticEnergy);

    rfc->setPhaseModelName(Attributes::getString(itsAttr[PHASE_MODEL]));
    rfc->setAmplitudeModelName(Attributes::getString(itsAttr[AMPLITUDE_MODEL]));
    rfc->setFrequencyModelName(Attributes::getString(itsAttr[FREQUENCY_MODEL]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(rfc);
}
