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
#ifndef OPAL_Quit_HH
#define OPAL_Quit_HH

#include "AbstractObjects/Action.h"


class Quit: public Action {

public:

    /// Exemplar constructor.
    Quit();

    virtual ~Quit();

    /// Make clone.
    virtual Quit* clone(const std::string& name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    Quit(const Quit&);
    void operator=(const Quit&);

    // Clone constructor.
    Quit(const std::string& name, Quit* parent);
};

#endif // OPAL_Quit_HH
