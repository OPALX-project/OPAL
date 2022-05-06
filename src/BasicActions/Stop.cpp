//
// Class Stop
//   The class for OPAL STOP command.
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
#include "BasicActions/Stop.h"


Stop::Stop(): Action(0, "STOP",
                         "The \"STOP\" statement terminates program execution "
                         "or reading of a called file.")
{}


Stop::Stop(const std::string& name, Stop* parent):
    Action(name, parent)
{}


Stop::~Stop()
{}


Stop* Stop::clone(const std::string& name) {
    return new Stop(name, this);
}


void Stop::execute()
{}
