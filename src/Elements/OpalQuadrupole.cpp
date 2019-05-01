// ------------------------------------------------------------------------
// $RCSfile: OpalQuadrupole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalQuadrupole
//   The class of OPAL Quadrupoles.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:07 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Elements/OpalQuadrupole.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MultipoleRep.h"
#include "ComponentWrappers/MultipoleWrapper.h"
#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include "Structure/ParticleMatterInteraction.h"

#include <cmath>
#include <iostream>
#include <sstream>


// Class OpalQuadrupole
// ------------------------------------------------------------------------

OpalQuadrupole::OpalQuadrupole():
    OpalElement(SIZE, "QUADRUPOLE",
                "The \"QUADRUPOLE\" element defines a Quadrupole."),
    parmatint_m(NULL) {
    itsAttr[K1] = Attributes::makeReal
                  ("K1", "Normalised upright quadrupole coefficient in m^(-2)");
    itsAttr[DK1] = Attributes::makeReal
                  ("DK1", "Normalised upright quadrupole coefficient error in m^(-2)");
    itsAttr[K1S] = Attributes::makeReal
                   ("K1S", "Normalised skew quadrupole coefficient in m^(-2)");

    itsAttr[DK1S] = Attributes::makeReal
                   ("DK1S", "Normalised skew quadrupole coefficient error in m^(-2)");
                   
    itsAttr[NSLICES] = Attributes::makeReal
                      ("NSLICES",
                      "The number of slices/ steps for this element in Map Tracking", 1);

    registerRealAttribute("K1");
    registerRealAttribute("DK1");
    registerRealAttribute("K1S");
    registerRealAttribute("DK1S");
    registerRealAttribute("NSLICES");

    registerOwnership();

    setElement((new MultipoleRep("QUADRUPOLE")));
}


OpalQuadrupole::OpalQuadrupole(const std::string &name, OpalQuadrupole *parent):
    OpalElement(name, parent),
    parmatint_m(NULL) {
    setElement((new MultipoleRep(name)));
}


OpalQuadrupole::~OpalQuadrupole() {
    if(parmatint_m)
        delete parmatint_m;
}


OpalQuadrupole *OpalQuadrupole::clone(const std::string &name) {
    return new OpalQuadrupole(name, this);
}


void OpalQuadrupole::print(std::ostream &os) const {
    OpalElement::print(os);
}


void OpalQuadrupole::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);

    // Get the desired field.
    const MultipoleWrapper *mult =
        dynamic_cast<const MultipoleWrapper *>(base.removeAlignWrapper());
    BMultipoleField field;

    // Get the desired field.
    if(flag == ERROR_FLAG) {
        field = mult->errorField();
    } else if(flag == ACTUAL_FLAG) {
        field = mult->getField();
    } else if(flag == IDEAL_FLAG) {
        field = mult->getDesign().getField();
    }

    double length = getLength();
    double scale = Physics::c / OpalData::getInstance()->getP0();
    if(length != 0.0) scale *= length;

    for(int order = 1; order <= field.order(); ++order) {
        std::ostringstream ss;
        ss << (order - 1) << std::ends;
        std::string orderString = ss.str();

        std::string normName = "K" + orderString + "L";
        registerRealAttribute(normName)->setReal(scale * field.normal(order));

        std::string skewName = "K" + orderString + "SL";
        registerRealAttribute(skewName)->setReal(scale * field.skew(order));

        scale *= double(order);
    }
}


void OpalQuadrupole::update() {
    OpalElement::update();

    MultipoleRep *quad =
        dynamic_cast<MultipoleRep *>(getElement()->removeWrappers());
    quad->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    double factor = OpalData::getInstance()->getP0() / Physics::c;

    BMultipoleField field;
    field.setNormalComponent(2, factor * Attributes::getReal(itsAttr[K1]));   // this is for the maps
    field.setSkewComponent(2, factor * Attributes::getReal(itsAttr[K1S]));    // this is for the maps
    quad->setField(field);
    quad->setNormalComponent(2, Attributes::getReal(itsAttr[K1]), Attributes::getReal(itsAttr[DK1]));
    quad->setSkewComponent(2, Attributes::getReal(itsAttr[K1S]), Attributes::getReal(itsAttr[DK1S]));
    quad->setNSlices(Attributes::getReal(itsAttr[NSLICES]));

    if(itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == NULL) {
        parmatint_m = (ParticleMatterInteraction::find(Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION])))->clone(getOpalName() + std::string("_parmatint"));
        parmatint_m->initParticleMatterInteractionHandler(*quad);
        quad->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(quad);
}
