//
// Class PSystem
//   The class for OPAL SYSTEM command in parallel environment.
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
#include "BasicActions/PSystem.h"

#include "Attributes/Attributes.h"
#include "Utility/IpplInfo.h"

#include <cstdlib>


PSystem::PSystem():
    Action(1, "PSYSTEM",
           "The \"PSYSTEM\" statement sends a command string to the "
           "operating system from all nodes.") {
    itsAttr[0] = Attributes::makeString
                 ("CMD", "A system command to be executed");

    registerOwnership(AttributeHandler::COMMAND);
}


PSystem::PSystem(const std::string& name, PSystem* parent):
    Action(name, parent)
{}


PSystem::~PSystem()
{}


PSystem* PSystem::clone(const std::string& name) {
    return new PSystem(name, this);
}


void PSystem::execute() {
    std::string command = Attributes::getString(itsAttr[0]);

    int res = system(command.c_str());
    if (res!=0)
        ERRORMSG("SYSTEM call failed" << endl);
}

void PSystem::parse(Statement& statement) {
    parseShortcut(statement);
}
