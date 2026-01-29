// ------------------------------------------------------------------------
// $RCSfile: Action.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Action
//   The base class for all OPAL action commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class Action
// ------------------------------------------------------------------------

Action::~Action()
{}


bool Action::canReplaceBy(Object *) {
    return true;
}


const std::string Action::getCategory() const {
    return "ACTION";
}


bool Action::shouldTrace() const {
    return true;

}

bool Action::shouldUpdate() const {
    return true;
}


Action::Action(const std::string &name, Action *parent):
    Object(name, parent)
{}


Action::Action(int size, const char *name, const char *help):
    Object(size, name, help)
{}
