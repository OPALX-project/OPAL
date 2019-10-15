// ------------------------------------------------------------------------
// $RCSfile: OpalBeamStripping.cpp,v $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalBeamStripping
//   The class of OPAL Cyclotron beam stripping.
//
// ------------------------------------------------------------------------
// $Date: 2018/11 $
// $Author: PedroCalvo$
// ------------------------------------------------------------------------

#include "Attributes/Attributes.h"
#include "BeamlineCore/BeamStrippingRep.h"
#include "Elements/OpalBeamStripping.h"
#include "Physics/Physics.h"
#include "Structure/ParticleMatterInteraction.h"

using Physics::pi;

// Class OpalBeamStripping
// ------------------------------------------------------------------------

OpalBeamStripping::OpalBeamStripping():
    OpalElement(SIZE, "BEAMSTRIPPING",
                "The \"BEAMSTRIPPING\" element defines a beam stripping interaction"),
    parmatint_m(NULL) {
    itsAttr[PRESSURE]        = Attributes::makeReal
        ("PRESSURE", " Pressure in the accelerator, [mbar]");
    itsAttr[TEMPERATURE]  = Attributes::makeReal
        ("TEMPERATURE", " Temperature of the accelerator, [K]");
    itsAttr[PMAPFN]   = Attributes::makeString
        ("PMAPFN", "Filename for the Pressure fieldmap");
    itsAttr[PSCALE]   = Attributes::makeReal
        ("PSCALE", "Scale factor for the P-field", 1.0);
    itsAttr[GAS]          = Attributes::makeString
        ("GAS", "The composition of residual gas");
    itsAttr[STOP]         = Attributes::makeBool
        ("STOP", "Option Whether stop tracking after beam stripping. Default: true", true);
    
    registerRealAttribute("PRESSURE");
    registerRealAttribute("TEMPERATURE");
    registerStringAttribute("PMAPFN");
    registerRealAttribute("PSCALE");
    registerStringAttribute("GAS");
    
    registerOwnership();
    
    setElement((new BeamStrippingRep("BEAMSTRIPPING"))->makeAlignWrapper());
}


OpalBeamStripping::OpalBeamStripping(const std::string &name, OpalBeamStripping *parent):
    OpalElement(name, parent),
    parmatint_m(NULL) {
    setElement((new BeamStrippingRep(name))->makeAlignWrapper());
}


OpalBeamStripping::~OpalBeamStripping() {
    delete parmatint_m;
}


OpalBeamStripping *OpalBeamStripping::clone(const std::string &name) {
    return new OpalBeamStripping(name, this);
}


void OpalBeamStripping::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
}


void OpalBeamStripping::update() {
    OpalElement::update();

    BeamStrippingRep *bstp =
        dynamic_cast<BeamStrippingRep *>(getElement()->removeWrappers());

    double pressure     = Attributes::getReal(itsAttr[PRESSURE]);
    double temperature  = Attributes::getReal(itsAttr[TEMPERATURE]);
    std::string pmap    = Attributes::getString(itsAttr[PMAPFN]);
    double pscale       = Attributes::getReal(itsAttr[PSCALE]);
    bool   stop         = Attributes::getBool(itsAttr[STOP]);
    std::string gas     = Attributes::getString(itsAttr[GAS]);

    bstp->setPressure(pressure);
    bstp->setTemperature(temperature);
    bstp->setPressureMapFN(pmap);
    bstp->setPScale(pscale);
    bstp->setStop(stop);
    bstp->setResidualGas(gas);

    if(itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == NULL) {
        parmatint_m = (ParticleMatterInteraction::find(Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION])))->clone(getOpalName() + std::string("_parmatint"));
        parmatint_m->initParticleMatterInteractionHandler(*bstp);
        bstp->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(bstp);
}