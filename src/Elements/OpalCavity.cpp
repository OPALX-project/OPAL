// ------------------------------------------------------------------------
// $RCSfile: OpalCavity.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalCavity
//   The class of OPAL RF cavities.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalCavity.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/RFCavityRep.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "Structure/OpalWake.h"
#include "Structure/BoundaryGeometry.h"
#include "Physics/Physics.h"

extern Inform *gmsg;

// Class OpalCavity
// ------------------------------------------------------------------------

OpalCavity::OpalCavity():
    OpalElement(SIZE, "RFCAVITY",
                "The \"RFCAVITY\" element defines an RF cavity."),
    owk_m(NULL),
    obgeo_m(NULL) {
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
    itsAttr[HARMON] = Attributes::makeReal
                      ("HARMON", "Harmonic number");
    itsAttr[BETARF] = Attributes::makeReal
                      ("BETRF", "beta_RF");
    itsAttr[PG] = Attributes::makeReal
                  ("PG", "RF power in MW");
    itsAttr[ZSHUNT] = Attributes::makeReal
                      ("SHUNT", "Shunt impedance in MOhm");
    itsAttr[TFILL] = Attributes::makeReal
                     ("TFILL", "Fill time in microseconds");
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

    registerRealAttribute("VOLT");
    registerRealAttribute("DVOLT");
    registerRealAttribute("FREQ");
    registerRealAttribute("LAG");
    registerRealAttribute("DLAG");
    registerStringAttribute("FMAPFN");
    registerStringAttribute("GEOMETRY");
    registerRealAttribute("RMIN");
    registerRealAttribute("RMAX");
    registerRealAttribute("ANGLE");
    registerRealAttribute("PDIS");
    registerRealAttribute("GAPWIDTH");
    registerRealAttribute("PHI0");
    registerRealAttribute("DESIGNENERGY");

    // attibutes for timedependent values
    registerStringAttribute("PHASE_MODEL");
    registerStringAttribute("AMPLITUDE_MODEL");
    registerStringAttribute("FREQUENCY_MODEL");

    registerOwnership();

    setElement((new RFCavityRep("RFCAVITY"))->makeAlignWrapper());
}


OpalCavity::OpalCavity(const std::string &name, OpalCavity *parent):
    OpalElement(name, parent),
    owk_m(NULL),
    obgeo_m(NULL) {
    setElement((new RFCavityRep(name))->makeAlignWrapper());
}


OpalCavity::~OpalCavity() {
    if(owk_m)
        delete owk_m;
}


OpalCavity *OpalCavity::clone(const std::string &name) {
    return new OpalCavity(name, this);
}


void OpalCavity::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    if(flag != ERROR_FLAG) {
        const RFCavityRep *rfc =
            dynamic_cast<const RFCavityRep *>(base.removeWrappers());
        attributeRegistry["VOLT"]->setReal(rfc->getAmplitude());
        attributeRegistry["FREQ"]->setReal(rfc->getFrequency());
        attributeRegistry["LAG"]->setReal(rfc->getPhase());
        attributeRegistry["FMAPFN"]->setString(rfc->getFieldMapFN());
    }
}


void OpalCavity::update() {
    OpalElement::update();

    using Physics::two_pi;
    RFCavityRep *rfc =
        dynamic_cast<RFCavityRep *>(getElement()->removeWrappers());

    double length = Attributes::getReal(itsAttr[LENGTH]);
    double peak  = Attributes::getReal(itsAttr[VOLT]);
    double peakError  = Attributes::getReal(itsAttr[DVOLT]);
    double phase  = Attributes::getReal(itsAttr[LAG]);
    double phaseError  = Attributes::getReal(itsAttr[DLAG]);
    double freq   = 1e6 * two_pi * Attributes::getReal(itsAttr[FREQ]);
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

    if(itsAttr[WAKEF] && owk_m == NULL) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + std::string("_wake"));
        owk_m->initWakefunction(*rfc);
        rfc->setWake(owk_m->wf_m);
    }

    if(itsAttr[GEOMETRY] && obgeo_m == NULL) {
        obgeo_m = (BoundaryGeometry::find(Attributes::getString(itsAttr[GEOMETRY])))->clone(getOpalName() + std::string("_geometry"));
        if(obgeo_m) {
            rfc->setBoundaryGeometry(obgeo_m);
        }
    }

    rfc->setElementLength(length);

    rfc->setAmplitude(1e6 * peak);
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
    rfc->setComponentType(type);
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