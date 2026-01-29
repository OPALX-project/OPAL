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
#ifndef OPAL_Call_HH
#define OPAL_Call_HH

#include "AbstractObjects/Action.h"


class Call: public Action {

public:

    /// Exemplar constructor.
    Call();

    virtual ~Call();

    /// Make clone.
    virtual Call* clone(const std::string& name);

    /// Execute the command.
    virtual void execute();

    /// Parse command (special for one-attribute command).
    virtual void parse(Statement&);

private:

    // Not implemented.
    Call(const Call&);
    void operator=(const Call&);

    // Clone constructor.
    Call(const std::string& name, Call* parent);
};

#endif // OPAL_Call_HH
