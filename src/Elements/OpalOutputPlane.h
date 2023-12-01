// Copyright (c) 2023, Chris Rogers
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

#ifndef OPAL_OpalOutputPlane_H
#define OPAL_OpalOutputPlane_H

#include "Elements/OpalElement.h"

// Class OpalOutputPlane
// ------------------------------------------------------------------------
/// Interface for output plane element.


class OpalOutputPlane: public OpalElement {

public:

    /// The attributes of class OpalOutputPlane.
    enum {
        CENTRE = COMMON, // plane centre position (3-vector)
        NORMAL,          // plane normal (3-vector)
        WIDTH,           // maximum allowed horizontal displacement from centre, in global x-y plane
        HEIGHT,          // maximum allowed vertical displacement from centre, in global vertical direction
        RADIUS,          // maximum allowed radial displacement from centre

        XSTART,            // "PROBE"-style placement - horizontal start position
        XEND,              // "PROBE"-style placement - horizontal end position
        YSTART,            // "PROBE"-style placement - vertical start position
        YEND,              // "PROBE"-style placement - vertical end position

        PLACEMENT_STYLE,   // "PROBE" or "CENTRE_NORMAL"

        ALGORITHM,         // algorithm used to find crossing-point
        TOLERANCE,         // tolerance on position estimate
        REFERENCE_ALIGNMENT_PARTICLE, // if true, centres on reference particle
        VERBOSE,           // set to 0 - 4 to get more verbose output

        SIZE
    };
    /// Exemplar constructor.
    OpalOutputPlane();

    virtual ~OpalOutputPlane();

    /// Make clone.
    virtual OpalOutputPlane *clone(const std::string &name);

    /// Fill in all registered attributes.
    //virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC septum.
    virtual void update();

private:

    // Not implemented.
    OpalOutputPlane(const OpalOutputPlane &);
    void operator=(const OpalOutputPlane &);

    // Clone constructor.
    OpalOutputPlane(const std::string &name, OpalOutputPlane *parent);
};

#endif // OPAL_OpalOutputPlane_H
