//
// Class Title
//   The class for OPAL TITLE command.
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
#include "BasicActions/Title.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"


Title::Title():
    Action(1, "TITLE",
           "The \"TITLE\" statement defines a new page title for subsequent "
           "output.") {
    itsAttr[0] = Attributes::makeString
                 ("STRING", "The title string");

    registerOwnership(AttributeHandler::STATEMENT);
}


Title::Title(const std::string& name, Title* parent):
    Action(name, parent)
{}


Title::~Title()
{}


Title* Title::clone(const std::string& name) {
    return new Title(name, this);
}


void Title::execute() {
    std::string title = Attributes::getString(itsAttr[0]);
    OpalData::getInstance()->storeTitle(title);
}


void Title::parse(Statement& statement) {
    parseShortcut(statement);
}
