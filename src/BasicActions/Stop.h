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
#ifndef OPAL_Stop_HH
#define OPAL_Stop_HH

#include "AbstractObjects/Action.h"


class Stop: public Action {

public:

    /// Exemplar constructor.
    Stop();

    virtual ~Stop();

    /// Make clone.
    virtual Stop* clone(const std::string& name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    Stop(const Stop&);
    void operator=(const Stop&);

    // Clone constructor.
    Stop(const std::string& name, Stop* parent);
};

#endif // OPAL_Stop_HH
