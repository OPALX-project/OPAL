// ------------------------------------------------------------------------
// $RCSfile: Definition.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Definition
//   The abstract base class for all OPAL definitions.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Definition.h"


// Class Definition
// ------------------------------------------------------------------------

Definition::~Definition()
{}


const std::string Definition::getCategory() const {
    return "DEFINITION";
}


bool Definition::shouldTrace() const {
    return false;

}

bool Definition::shouldUpdate() const {
    return false;
}


Definition::Definition(int size, const char *name, const char *help):
    Object(size, name, help)
{}


Definition::Definition(const std::string &name, Definition *parent):
    Object(name, parent)
{}
