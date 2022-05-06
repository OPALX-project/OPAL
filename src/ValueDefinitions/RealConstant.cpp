//
// Class RealConstant
//   The REAL CONSTANT definition.
//
// Copyright (c) 2000 - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "ValueDefinitions/RealConstant.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "OPALconfig.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Utilities/Options.h"
#include "Utility/IpplInfo.h"

#include <cmath>
#include <iostream>


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
    opal->create(new RealConstant("RADDEG", this, Units::rad2deg));
    opal->create(new RealConstant("DEGRAD", this, Units::deg2rad));
    opal->create(new RealConstant("E",      this, Physics::e));

    opal->create(new RealConstant("MUMASS",    this, Physics::m_mu));
    opal->create(new RealConstant("EMASS",     this, Physics::m_e));
    opal->create(new RealConstant("PMASS",     this, Physics::m_p));
    opal->create(new RealConstant("HMMASS",    this, Physics::m_hm));
    opal->create(new RealConstant("H2PMASS",   this, Physics::m_h2p));
    opal->create(new RealConstant("DMASS",     this, Physics::m_d));
    opal->create(new RealConstant("ALPHAMASS", this, Physics::m_alpha));
    opal->create(new RealConstant("CMASS",     this, Physics::m_c));
    opal->create(new RealConstant("XEMASS",    this, Physics::m_xe));
    opal->create(new RealConstant("UMASS",     this, Physics::m_u));

    opal->create(new RealConstant("CLIGHT", this, Physics::c));
    opal->create(new RealConstant("AMU",    this, Physics::amu));

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