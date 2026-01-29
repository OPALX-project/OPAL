//
// Class System
//   The class for OPAL SYSTEM command on a single core.
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
#include "BasicActions/System.h"

#include "Attributes/Attributes.h"
#include "Utility/IpplInfo.h"

#include <cstdlib>

System::System():
    Action(1, "SYSTEM",
           "The \"SYSTEM\" statement sends a command string to the "
           "operating system.") {
    itsAttr[0] = Attributes::makeString
                 ("CMD", "A system command to be executed");

    registerOwnership(AttributeHandler::COMMAND);
}


System::System(const std::string& name, System* parent):
    Action(name, parent)
{}


System::~System()
{}


System* System::clone(const std::string& name) {
    return new System(name, this);
}


void System::execute() {
    if (Ippl::myNode() == 0) {
        std::string command = Attributes::getString(itsAttr[0]);

        int res = system(command.c_str());
        if (res!=0)
            ERRORMSG("SYSTEM call failed" << endl);
    }
}

void System::parse(Statement& statement) {
    parseShortcut(statement);
}
