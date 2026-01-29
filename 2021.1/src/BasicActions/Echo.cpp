//
// Class Echo
//   The class for OPAL ECHO command.
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
#include "BasicActions/Echo.h"

#include "Attributes/Attributes.h"
#include "Utility/IpplInfo.h"

#include <iostream>


Echo::Echo():
    Action(1, "ECHO",
           "The \"ECHO\" statement sends a message to the ECHO file.") {
    itsAttr[0] = Attributes::makeString("MESSAGE", "The message to be sent.");

    registerOwnership(AttributeHandler::STATEMENT);
}


Echo::Echo(const std::string& name, Echo* parent):
    Action(name, parent)
{}


Echo::~Echo()
{}


Echo* Echo::clone(const std::string& name) {
    return new Echo(name, this);
}


void Echo::execute() {
    if (Ippl::myNode() == 0)
        std::cerr << Attributes::getString(itsAttr[0]) << std::endl;
}


void Echo::parse(Statement& statement) {
    parseShortcut(statement);
}
