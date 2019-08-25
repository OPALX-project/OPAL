// ------------------------------------------------------------------------
// $RCSfile: StringConstant.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: StringConstant
//   Implements a OPAL STRING_CONSTANT definition.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:26:42 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "ValueDefinitions/StringConstant.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Utilities/Util.h"
#include <iostream>


// Class StringConstant
// ------------------------------------------------------------------------

StringConstant::StringConstant():
    ValueDefinition(1, "STRING_CONSTANT",
                    "The \"STRING CONSTANT\" statement defines a global "
                    "string constant:\n"
                    "\tSTRING CONSTANT <name> = <String-expression>;\n") {
    itsAttr[0] = Attributes::makeString("VALUE", "The constant value");

    registerOwnership(AttributeHandler::STATEMENT);

    OpalData *opal = OpalData::getInstance();
    opal->create(new StringConstant("GITREVISION", this, Util::getGitRevision()));
}


StringConstant::StringConstant(const std::string &name, StringConstant *parent):
    ValueDefinition(name, parent)
{}


StringConstant::StringConstant(const std::string &name, StringConstant *parent, const std::string &value):
    ValueDefinition(name, parent)
{
    Attributes::setString(itsAttr[0], value);
    itsAttr[0].setReadOnly(true);
    builtin = true;
}


StringConstant::~StringConstant()
{}


bool StringConstant::canReplaceBy(Object *) {
    return false;
}


StringConstant *StringConstant::clone(const std::string &name) {
    return new StringConstant(name, this);
}



std::string StringConstant::getString() const {
    return Attributes::getString(itsAttr[0]);
}



void StringConstant::print(std::ostream &os) const {
    os << "STRING " << getOpalName() << '=' << itsAttr[0] << ';';
    os << std::endl;
}

void StringConstant::printValue(std::ostream &os) const {
    os << itsAttr[0];
}