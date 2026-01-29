//
// Class OpalProbe
//   The Probe element.
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
#ifndef OPAL_OpalProbe_HH
#define OPAL_OpalProbe_HH

#include "Elements/OpalElement.h"

class OpalWake;

class OpalProbe: public OpalElement {

public:

    /// The attributes of class OpalProbe.
    enum {
        XSTART = COMMON, // Start of x coordinate
        XEND,            // End of x coordinate
        YSTART,          // Start of y coordinate
        YEND,            // End of y coordinate
        WIDTH,           // Width of the probe
        STEP,            // Step size of the probe
        SIZE
    };
    /// Exemplar constructor.
    OpalProbe();

    virtual ~OpalProbe();

    /// Make clone.
    virtual OpalProbe *clone(const std::string &name);

    /// Update the embedded CLASSIC probe.
    virtual void update();

private:

    // Not implemented.
    OpalProbe(const OpalProbe &);
    void operator=(const OpalProbe &);

    // Clone constructor.
    OpalProbe(const std::string &name, OpalProbe *parent);

    OpalWake *owk_m;
};

#endif // OPAL_OpalProbe_HH
