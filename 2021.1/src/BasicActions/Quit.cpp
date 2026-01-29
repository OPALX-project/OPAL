//
// Class Quit
//   The class for OPAL QUIT command.
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
#include "BasicActions/Quit.h"


Quit::Quit(): Action(0, "QUIT",
                         "The \"QUIT\" statement terminates program execution "
                         "or reading of a called file.")
{}


Quit::Quit(const std::string& name, Quit* parent):
    Action(name, parent)
{}


Quit::~Quit()
{}


Quit* Quit::clone(const std::string& name) {
    return new Quit(name, this);
}


void Quit::execute()
{}
