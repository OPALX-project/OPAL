//
// Class OpalHKicker
//   The HKICKER element.
//   Note the sign convention:  A positive kick bend particles to positive x.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_OpalHKicker_HH
#define OPAL_OpalHKicker_HH

#include "Elements/OpalElement.h"

class OpalHKicker: public OpalElement {

public:

    /// The attributes of class OpalHKicker.
    enum {
        KICK = COMMON,  // The kicker strength.
        DESIGNENERGY,   // The mean kinetic energy at exit
        K0,             // The magnetic field
        SIZE
    };

    /// Exemplar constructor.
    OpalHKicker();

    virtual ~OpalHKicker();

    /// Make clone.
    virtual OpalHKicker *clone(const std::string &name);

    /// Update the embedded CLASSIC corrector.
    virtual void update();

private:

    // Not implemented.
    OpalHKicker(const OpalHKicker &);
    void operator=(const OpalHKicker &);

    // Clone constructor.
    OpalHKicker(const std::string &name, OpalHKicker *parent);
};

#endif // OPAL_OpalHKicker_HH
