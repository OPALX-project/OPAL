// ------------------------------------------------------------------------
// $RCSfile: OpalSextupole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalSextupole
//   The class of OPAL Sextupoles.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Elements/OpalSextupole.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MultipoleRep.h"
#include "ComponentWrappers/MultipoleWrapper.h"
#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include <cmath>
#include <iostream>
#include <sstream>


// Class OpalSextupole
// ------------------------------------------------------------------------

OpalSextupole::OpalSextupole():
    OpalElement(SIZE, "SEXTUPOLE",
                "The \"SEXTUPOLE\" element defines a Sextupole.") {
    itsAttr[K2] = Attributes::makeReal
                  ("K2", "Normalised upright sextupole coefficient in m^(-3)");
    itsAttr[DK2] = Attributes::makeReal
                  ("DK2", "Normalised upright sextupole coefficient error in m^(-3)");
    itsAttr[K2S] = Attributes::makeReal
                   ("K2S", "Normalised skew sextupole coefficient in m^(-3)");
    itsAttr[DK2S] = Attributes::makeReal
                   ("DK2S", "Normalised skew sextupole coefficient error in m^(-3)");

    registerOwnership();

    setElement((new MultipoleRep("SEXTUPOLE"))->makeWrappers());
}


OpalSextupole::OpalSextupole(const std::string &name, OpalSextupole *parent):
    OpalElement(name, parent) {
    setElement((new MultipoleRep(name))->makeWrappers());
}


OpalSextupole::~OpalSextupole()
{}


OpalSextupole *OpalSextupole::clone(const std::string &name) {
    return new OpalSextupole(name, this);
}



void OpalSextupole::print(std::ostream &os) const {
    OpalElement::print(os);
}


void OpalSextupole::
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


void OpalSextupole::update() {
    OpalElement::update();

    MultipoleRep *sext =
        dynamic_cast<MultipoleRep *>(getElement()->removeWrappers());
    sext->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    double factor = OpalData::getInstance()->getP0() / (Physics::c * 2.0);
    BMultipoleField field;
    field.setNormalComponent(3, factor * Attributes::getReal(itsAttr[K2]));
    field.setSkewComponent(3, factor * Attributes::getReal(itsAttr[K2S]));
    sext->setField(field);
    sext->setNormalComponent(3, Attributes::getReal(itsAttr[K2]), Attributes::getReal(itsAttr[DK2]));
    sext->setSkewComponent(3, Attributes::getReal(itsAttr[K2S]), Attributes::getReal(itsAttr[DK2S]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(sext);
}