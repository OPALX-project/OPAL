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
#ifndef OPAL_System_HH
#define OPAL_System_HH

#include "AbstractObjects/Action.h"


class System: public Action {

public:

    /// Exemplar constructor.
    System();

    virtual ~System();

    /// Make clone.
    virtual System* clone(const std::string& name);

    /// Execute the command.
    virtual void execute();

    /// Parse command (special for one-attribute command).
    virtual void parse(Statement&);

private:

    // Not implemented.
    System(const System&);
    void operator=(const System&);

    // Clone constructor.
    System(const std::string& name, System* parent);
};

#endif // OPAL_System_HH
