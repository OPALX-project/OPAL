//
// Class Call
//   The class for OPAL CALL command.
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
#include "BasicActions/Call.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "OpalParser/OpalParser.h"
#include "Parser/FileStream.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utility/IpplInfo.h"

#include <iostream>

extern Inform* gmsg;

Call::Call():
    Action(1, "CALL",
           "The \"CALL\" statement switches input temporarily to the "
           "named file.") {
    itsAttr[0] = Attributes::makeString
                 ("FILE", "Name of file to be read", "CALL");

    registerOwnership(AttributeHandler::STATEMENT);
}


Call::Call(const std::string& name, Call* parent):
    Action(name, parent)
{}


Call::~Call()
{}


Call* Call::clone(const std::string& name) {
    return new Call(name, this);
}


void Call::execute() {
    std::string file = Attributes::getString(itsAttr[0]);

    if (Options::info && Ippl::myNode() == 0) {
        *gmsg << "* Reading input stream '" << file
              << "' from \"CALL\" command.\n" << endl;
    }

    OpalParser().run(new FileStream(file));
}


void Call::parse(Statement& statement) {
    parseShortcut(statement);
}
