// ------------------------------------------------------------------------
// $RCSfile: OpalTravelingWave.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalTravelingWave
//   The class of OPAL RF cavities.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalTravelingWave.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/TravelingWaveRep.h"
#include "Structure/OpalWake.h"
#include "Physics/Physics.h"


// Class OpalTravelingWave
// ------------------------------------------------------------------------

OpalTravelingWave::OpalTravelingWave():
    OpalElement(SIZE, "TRAVELINGWAVE",
                "The \"TRAVELINGWAVE\" element defines a traveling wave structure."),
    owk_m(NULL) {
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
                      ("FMAPFN", "Filename for the fieldmap");
    itsAttr[FAST] = Attributes::makeBool
                    ("FAST", "Faster but less accurate", true);
    itsAttr[APVETO] = Attributes::makeBool
                    ("APVETO", "Do not use this cavity in the Autophase procedure", false);
    itsAttr[CAVITYTYPE] = Attributes::makeString
                          ("CAVITYTYPE", "STANDING or TRAVELING wave cavity in photoinjector and LINAC; SINGLEGAP or DOUBLEGAP cavity in cyclotron");
    itsAttr[NUMCELLS] = Attributes::makeReal
                        ("NUMCELLS", "Number of cells in a TW structure");
    itsAttr[DESIGNENERGY] = Attributes::makeReal
                            ("DESIGNENERGY", "the mean energy of the particles at exit", -1.0);
    itsAttr[MODE] = Attributes::makeReal
                     ("MODE", "The phase shift between neighboring cells in 2*pi", 1.0/3.0);

    registerRealAttribute("VOLT");
    registerRealAttribute("DVOLT");
    registerRealAttribute("FREQ");
    registerRealAttribute("LAG");
    registerRealAttribute("DLAG");
    registerStringAttribute("FMAPFN");
    registerStringAttribute("CAVITYTYPE");
    registerRealAttribute("NUMCELLS");
    registerRealAttribute("DESIGNENERGY");
    registerRealAttribute("MODE");

    registerOwnership();

    setElement((new TravelingWaveRep("TRAVELINGWAVE"))->makeAlignWrapper());
}


OpalTravelingWave::OpalTravelingWave(const std::string &name, OpalTravelingWave *parent):
    OpalElement(name, parent),
    owk_m(NULL) {
    setElement((new TravelingWaveRep(name))->makeAlignWrapper());
}


OpalTravelingWave::~OpalTravelingWave() {
    if(owk_m)
        delete owk_m;
}


OpalTravelingWave *OpalTravelingWave::clone(const std::string &name) {
    return new OpalTravelingWave(name, this);
}


void OpalTravelingWave::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    if(flag != ERROR_FLAG) {
        const TravelingWaveRep *rfc =
            dynamic_cast<const TravelingWaveRep *>(base.removeWrappers());
        attributeRegistry["VOLT"]->setReal(rfc->getAmplitude());
        attributeRegistry["DVOLT"]->setReal(rfc->getAmplitudeError());
        attributeRegistry["FREQ"]->setReal(rfc->getFrequency());
        attributeRegistry["LAG"]->setReal(rfc->getPhase());
        attributeRegistry["DLAG"]->setReal(rfc->getPhaseError());
        attributeRegistry["FMAPFN"]->setString(rfc->getFieldMapFN());
    }
}


void OpalTravelingWave::update() {
    OpalElement::update();

    TravelingWaveRep *rfc =
        dynamic_cast<TravelingWaveRep *>(getElement()->removeWrappers());

    double length = Attributes::getReal(itsAttr[LENGTH]);
    double vPeak  = Attributes::getReal(itsAttr[VOLT]);
    double vPeakError  = Attributes::getReal(itsAttr[DVOLT]);
    double phase  = Attributes::getReal(itsAttr[LAG]);
    double phaseError  = Attributes::getReal(itsAttr[DLAG]);
    double freq   = (1.0e6 * Physics::two_pi) * Attributes::getReal(itsAttr[FREQ]);
    std::string fmapfm = Attributes::getString(itsAttr[FMAPFN]);
    bool fast = Attributes::getBool(itsAttr[FAST]);
    bool apVeto = Attributes::getBool(itsAttr[APVETO]);

    //    std::string type = Attributes::getString(itsAttr[TYPE]);
    double kineticEnergy = Attributes::getReal(itsAttr[DESIGNENERGY]);

    rfc->setElementLength(length);
    rfc->setAmplitude(1.0e6 * vPeak);
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

    if(itsAttr[WAKEF] && owk_m == NULL) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + std::string("_wake"));
        owk_m->initWakefunction(*rfc);
        rfc->setWake(owk_m->wf_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(rfc);
}