//
// Class Select
//   The class for OPAL SELECT command.
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
#ifndef OPAL_Select_HH
#define OPAL_Select_HH 1

#include "AbstractObjects/Action.h"

class Beamline;


class Select: public Action {

public:

    /// Exemplar constructor.
    Select();

    virtual ~Select();

    /// Make clone.
    virtual Select* clone(const std::string& name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    Select(const Select&);
    void operator=(const Select&);

    // Clone constructor.
    Select(const std::string& name, Select* parent);

    // Do the selection.
    void select(const Beamline&);
};

#endif // OPAL_Select_H
