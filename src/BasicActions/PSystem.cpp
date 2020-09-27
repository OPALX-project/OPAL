// ------------------------------------------------------------------------
// $RCSfile: PSystem.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: PSystem
//   The class for the OPAL SYSTEM command.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/01/17 22:18:36 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "BasicActions/PSystem.h"

#include "Attributes/Attributes.h"

#include "Utility/IpplInfo.h"

#include <cstdlib>
// Class PSystem
// ------------------------------------------------------------------------

PSystem::PSystem():
    Action(1, "PSYSTEM",
           "The \"PSYSTEM\" statement sends a command string to the "
           "operating system from all nodes.") {
    itsAttr[0] = Attributes::makeString
                 ("CMD", "A system command to be executed");

    registerOwnership(AttributeHandler::COMMAND);
}


PSystem::PSystem(const std::string &name, PSystem *parent):
    Action(name, parent)
{}


PSystem::~PSystem()
{}


PSystem *PSystem::clone(const std::string &name) {
    return new PSystem(name, this);
}


void PSystem::execute() {
    std::string command = Attributes::getString(itsAttr[0]);

    int res = system(command.c_str());
    if (res!=0)
        ERRORMSG("SYSTEM call failed" << endl);
}

void PSystem::parse(Statement &statement) {
    parseShortcut(statement);
}