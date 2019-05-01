// ------------------------------------------------------------------------
// $RCSfile: OpalOctupole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalOctupole
//   The class of OPAL Octupoles.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:07 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Elements/OpalOctupole.h"
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

// Class OpalOctupole
// ------------------------------------------------------------------------

OpalOctupole::OpalOctupole():
    OpalElement(SIZE, "OCTUPOLE",
                "The \"OCTUPOLE\" element defines a Octupole.") {
    itsAttr[K3] = Attributes::makeReal
                  ("K3", "Normalised upright octupole coefficient in m^(-4)");
    itsAttr[DK3] = Attributes::makeReal
                  ("DK3", "Normalised upright octupole coefficient error in m^(-4)");
    itsAttr[K3S] = Attributes::makeReal
                   ("K3S", "Normalised skew octupole coefficient in m^(-4)");
    itsAttr[DK3S] = Attributes::makeReal
                   ("DK3S", "Normalised skew octupole coefficient error in m^(-4)");

    registerOwnership();

    setElement((new MultipoleRep("OCTUPOLE"))->makeWrappers());
}


OpalOctupole::OpalOctupole(const std::string &name, OpalOctupole *parent):
    OpalElement(name, parent) {
    setElement((new MultipoleRep(name))->makeWrappers());
}


OpalOctupole::~OpalOctupole()
{}


OpalOctupole *OpalOctupole::clone(const std::string &name) {
    return new OpalOctupole(name, this);
}


void OpalOctupole::print(std::ostream &os) const {
    OpalElement::print(os);
}


void OpalOctupole::
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


void OpalOctupole::update() {
    OpalElement::update();

    MultipoleRep *oct =
        dynamic_cast<MultipoleRep *>(getElement()->removeWrappers());
    oct->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    double factor = OpalData::getInstance()->getP0() / (Physics::c * 6.0);
    BMultipoleField field;
    field.setNormalComponent(4, factor * Attributes::getReal(itsAttr[K3]));
    field.setSkewComponent(4, factor * Attributes::getReal(itsAttr[K3S]));
    oct->setField(field);

    oct->setNormalComponent(4, Attributes::getReal(itsAttr[K3]), Attributes::getReal(itsAttr[DK3]));
    oct->setSkewComponent(4, Attributes::getReal(itsAttr[K3S]), Attributes::getReal(itsAttr[DK3S]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(oct);
}