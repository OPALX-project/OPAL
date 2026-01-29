//
// Class SampleCmd
//   This class defines the SAMPLE command.
//
// Copyright (c) 2018, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
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
#ifndef OPAL_SampleCmd_HH
#define OPAL_SampleCmd_HH

#include "AbstractObjects/Action.h"

#include <string>

class SampleCmd: public Action {

public:

    /// Exemplar constructor.
    SampleCmd();

    virtual ~SampleCmd();

    /// Make clone.
    virtual SampleCmd *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    SampleCmd(const SampleCmd &)      = delete;
    void operator=(const SampleCmd &) = delete;

    // Clone constructor.
    SampleCmd(const std::string &name, SampleCmd *parent);

    void stashEnvironment();
    void popEnvironment();
};

#endif