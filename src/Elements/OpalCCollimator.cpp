// ------------------------------------------------------------------------
// $RCSfile: OpalCCollimator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalCCollimator
//   The class of OPAL Cyclotron collimators.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann, Jianjun Yang $
//
// ------------------------------------------------------------------------

#include "Elements/OpalCCollimator.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/CCollimatorRep.h"
#include "Structure/ParticleMatterInteraction.h"
#include "Physics/Physics.h"


// Class OpalCCollimator
// ------------------------------------------------------------------------

OpalCCollimator::OpalCCollimator():
    OpalElement(SIZE, "CCOLLIMATOR",
                "The \"CCOLLIMATOR\" element defines a rectangular-shape cyclotron collimator"),
    parmatint_m(NULL) {
    itsAttr[XSTART] = Attributes::makeReal
                      ("XSTART", " Start of x coordinate [mm]");
    itsAttr[XEND] = Attributes::makeReal
                    ("XEND", " End of x coordinate, [mm]");
    itsAttr[YSTART] = Attributes::makeReal
                      ("YSTART", "Start of y coordinate, [mm]");
    itsAttr[YEND] = Attributes::makeReal
                    ("YEND", "End of y coordinate, [mm]");
    itsAttr[ZSTART] = Attributes::makeReal
      ("ZSTART", "Start of vertical coordinate, [mm], default value: -100",-100.0);
    itsAttr[ZEND] = Attributes::makeReal
      ("ZEND", "End of vertical coordinate, [mm], default value: 100", 100.0);
    itsAttr[WIDTH] = Attributes::makeReal
                     ("WIDTH", "Width of the collimator [mm]");
    itsAttr[OUTFN] = Attributes::makeString
                     ("OUTFN", "Output filename");

    registerRealAttribute("XSTART");
    registerRealAttribute("XEND");
    registerRealAttribute("YSTART");
    registerRealAttribute("YEND");
    registerRealAttribute("ZSTART");
    registerRealAttribute("ZEND");
    registerRealAttribute("WIDTH");
    registerStringAttribute("OUTFN");

    registerOwnership();

    setElement((new CCollimatorRep("CCOLLIMATOR"))->makeAlignWrapper());
}


OpalCCollimator::OpalCCollimator(const std::string &name, OpalCCollimator *parent):
    OpalElement(name, parent),
    parmatint_m(NULL) {
    setElement((new CCollimatorRep(name))->makeAlignWrapper());
}


OpalCCollimator::~OpalCCollimator() {
    if(parmatint_m)
        delete parmatint_m;
}


OpalCCollimator *OpalCCollimator::clone(const std::string &name) {
    return new OpalCCollimator(name, this);
}


void OpalCCollimator::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
}


void OpalCCollimator::update() {
    OpalElement::update();

    CCollimatorRep *coll =
        dynamic_cast<CCollimatorRep *>(getElement()->removeWrappers());
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double xstart = Attributes::getReal(itsAttr[XSTART]);
    double xend = Attributes::getReal(itsAttr[XEND]);
    double ystart = Attributes::getReal(itsAttr[YSTART]);
    double yend = Attributes::getReal(itsAttr[YEND]);
    double zstart = Attributes::getReal(itsAttr[ZSTART]);
    double zend = Attributes::getReal(itsAttr[ZEND]);
    double width = Attributes::getReal(itsAttr[WIDTH]);
    coll->setElementLength(length);
    coll->setDimensions(xstart, xend, ystart, yend, zstart, zend, width);
    coll->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    if(itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == NULL) {
        parmatint_m = (ParticleMatterInteraction::find(Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION])))->clone(getOpalName() + std::string("_parmatint"));
        parmatint_m->initParticleMatterInteractionHandler(*coll);
        coll->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(coll);
}