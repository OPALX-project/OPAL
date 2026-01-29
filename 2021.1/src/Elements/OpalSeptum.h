//
// Class OpalSeptum
//   The Septum element.
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
#ifndef OPAL_OpalSeptum_HH
#define OPAL_OpalSeptum_HH

#include "Elements/OpalElement.h"

class OpalWake;

class OpalSeptum: public OpalElement {

public:

    /// The attributes of class OpalSeptum.
    enum {
        XSTART = COMMON, // Start of x coordinate
        XEND,            // End of x coordinate
        YSTART,          // Start of y coordinate
        YEND,            // End of y coordinate
        WIDTH,           // Width of the septum
        SIZE
    };
    /// Exemplar constructor.
    OpalSeptum();

    virtual ~OpalSeptum();

    /// Make clone.
    virtual OpalSeptum *clone(const std::string &name);

    /// Update the embedded CLASSIC septum.
    virtual void update();

private:

    // Not implemented.
    OpalSeptum(const OpalSeptum &);
    void operator=(const OpalSeptum &);

    // Clone constructor.
    OpalSeptum(const std::string &name, OpalSeptum *parent);

    OpalWake *owk_m;
};

#endif // OPAL_OpalSeptum_HH
