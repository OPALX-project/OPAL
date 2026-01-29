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
#ifndef OPAL_PSystem_HH
#define OPAL_PSystem_HH

#include "AbstractObjects/Action.h"


class PSystem: public Action {

public:

    /// Exemplar constructor.
    PSystem();

    virtual ~PSystem();

    /// Make clone.
    virtual PSystem* clone(const std::string& name);

    /// Execute the command.
    virtual void execute();

    /// Parse command (special for one-attribute command).
    virtual void parse(Statement&);

private:

    // Not implemented.
    PSystem(const PSystem&);
    void operator=(const PSystem&);

    // Clone constructor.
    PSystem(const std::string& name, PSystem* parent);
};

#endif // OPAL_PSystem_HH
