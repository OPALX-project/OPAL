// ------------------------------------------------------------------------
// $RCSfile: ValueDefinition.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ValueDefinition
//   The abstract base class for all OPAL value definitions.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:35 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/ValueDefinition.h"
#include "AbstractObjects/Object.h"
#include "Parser/Statement.h"
#include "Utilities/OpalException.h"


// Class ValueDefinition
// ------------------------------------------------------------------------

ValueDefinition::~ValueDefinition()
{}


const std::string ValueDefinition::getCategory() const {
    return "VARIABLE";
}


bool ValueDefinition::shouldTrace() const {
    return false;

}

bool ValueDefinition::shouldUpdate() const {
    return false;
}


bool ValueDefinition::getBool() const {
    throw OpalException("ValueDefinition::getBool()",
                        "You cannot get a logical value from \"" +
                        getOpalName() + "\".");
}


bool ValueDefinition::getBoolComponent(int) const {
    throw OpalException("ValueDefinition::getBoolComponent()",
                        "You cannot get an indexed logical value from \"" +
                        getOpalName() + "\".");
}


double ValueDefinition::getReal() const {
    throw OpalException("ValueDefinition::getReal()",
                        "You cannot get a real value from \"" +
                        getOpalName() + "\".");
}


double ValueDefinition::getRealComponent(int) const {
    throw OpalException("ValueDefinition::getReal()",
                        "You cannot get an indexed real value from \"" +
                        getOpalName() + "\".");
}


std::string ValueDefinition::getString() const {
    throw OpalException("ValueDefinition::getString()",
                        "You cannot get a string value from \"" +
                        getOpalName() + "\".");
}


std::string ValueDefinition::getStringComponent(int) const {
    throw OpalException("ValueDefinition::getString()",
                        "You cannot get an indexed string value from \"" +
                        getOpalName() + "\".");
}


Attribute &ValueDefinition::value() {
    return itsAttr[0];
}


const Attribute &ValueDefinition::value() const {
    return itsAttr[0];
}


ValueDefinition::ValueDefinition
(int size, const char *name, const char *help):
    Object(size, name, help)
{}


ValueDefinition::ValueDefinition(const std::string &name, ValueDefinition *parent):
    Object(name, parent)
{}
