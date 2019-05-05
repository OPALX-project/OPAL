// ------------------------------------------------------------------------
// $RCSfile: RealConstant.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RealConstant
//   Implements a REAL_CONSTANT definition.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:26:42 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "ValueDefinitions/RealConstant.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Utilities/Options.h"
#include "Physics/Physics.h"
#include "OPALconfig.h"

#include "Utility/IpplInfo.h"

#include <cmath>
#include <iostream>

// Class RealConstant
// ------------------------------------------------------------------------

RealConstant::RealConstant():
    ValueDefinition(1, "REAL_CONSTANT",
                    "The \"REAL CONSTANT\" statement defines a global "
                    "real constant:\n"
                    "\tREAL CONSTANT <name> = <real-expression>;\n") {
    itsAttr[0] = Attributes::makeReal("VALUE", "The constant value", 0.0);

    registerOwnership(AttributeHandler::STATEMENT);

    // Define the standard constants.
    OpalData *opal = OpalData::getInstance();
    opal->create(new RealConstant("PI",     this, Physics::pi));
    opal->create(new RealConstant("TWOPI",  this, Physics::two_pi));
    opal->create(new RealConstant("RADDEG", this, 180.0 / Physics::pi));
    opal->create(new RealConstant("DEGRAD", this, Physics::pi / 180.0));
    opal->create(new RealConstant("E",      this, Physics::e));

    opal->create(new RealConstant("EMASS",  this, Physics::m_e));
    opal->create(new RealConstant("PMASS",  this, Physics::m_p));
    opal->create(new RealConstant("HMMASS", this, Physics::m_hm));
    opal->create(new RealConstant("UMASS", this, Physics::m_u));
    opal->create(new RealConstant("CMASS", this, Physics::m_c));
    opal->create(new RealConstant("MMASS", this, Physics::m_mu));
    opal->create(new RealConstant("DMASS", this, Physics::m_d));
    opal->create(new RealConstant("XEMASS", this, Physics::m_xe));

    opal->create(new RealConstant("CLIGHT", this, Physics::c));

    opal->create(new RealConstant("OPALVERSION", this, OPAL_VERSION_MAJOR * 10000
                                  + OPAL_VERSION_MINOR * 100
                                  + OPAL_VERSION_PATCH));
    opal->create(new RealConstant("RANK", this, Ippl::myNode()));
}


RealConstant::RealConstant(const std::string &name, RealConstant *parent):
    ValueDefinition(name, parent)
{}


RealConstant::RealConstant(const std::string &name, RealConstant *parent,
                           double value):
    ValueDefinition(name, parent) {
    Attributes::setReal(itsAttr[0], value);
    itsAttr[0].setReadOnly(true);
    builtin = true;
}


RealConstant::~RealConstant()
{}


bool RealConstant::canReplaceBy(Object *) {
    return false;
}


RealConstant *RealConstant::clone(const std::string &name) {
    return new RealConstant(name, this);
}


double RealConstant::getReal() const {
    return Attributes::getReal(itsAttr[0]);
}


void RealConstant::print(std::ostream &os) const {
    os << "REAL CONST " << getOpalName() << '=' << itsAttr[0] << ';';
    os << std::endl;
}

void RealConstant::printValue(std::ostream &os) const {
    os << itsAttr[0];
}