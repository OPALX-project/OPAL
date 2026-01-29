//
// Class Help
//   The class for OPAL HELP command.
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
#include "BasicActions/Help.h"

#include "AbstractObjects/Object.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Utility/IpplInfo.h"

#include <iostream>


Help::Help():
    Action(1, "HELP",
           "The \"HELP\" statement displays the purpose and attribute "
           "types of an object.") {
    itsAttr[0] = Attributes::makeString
                 ("NAME", "Name of object for which help is wanted");

    registerOwnership(AttributeHandler::STATEMENT);
}


Help::Help(const std::string& name, Help* parent):
    Action(name, parent)
{}


Help::~Help()
{}


Help* Help::clone(const std::string& name) {
    return new Help(name, this);
}


void Help::execute() {
    if (itsAttr[0]) {
        std::string name = Attributes::getString(itsAttr[0]);

        if (Object *object = OpalData::getInstance()->find(name)) {
            object->printHelp(std::cerr);
        } else {
            ERRORMSG("\n" << *this << "Unknown object \"" << name << "\".\n"
                     << endl) ;
        }
    } else {
        printHelp(std::cerr);
    }
}


void Help::parse(Statement& statement) {
    parseShortcut(statement);
}
