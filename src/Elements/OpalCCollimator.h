//
// Class OpalCCollimator
//   The CCOLLIMATOR element.
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
#ifndef OPAL_OpalCCollimator_HH
#define OPAL_OpalCCollimator_HH

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;

class OpalCCollimator: public OpalElement {

public:

    /// The attributes of class OpalCCollimator.
    enum {
        XSTART = COMMON,  // Start of x coordinate
        XEND,             // End of x coordinate
        YSTART,           // Start of y coordinate
        YEND,             // End of y coordinate
        ZSTART,           // Top boundary
        ZEND,             // Bottom boundary
        WIDTH,            // The width of collimator
        SIZE
    };

    /// Exemplar constructor.
    OpalCCollimator();

    virtual ~OpalCCollimator();

    /// Make clone.
    virtual OpalCCollimator* clone(const std::string& name);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalCCollimator(const OpalCCollimator&);
    void operator=(const OpalCCollimator&);

    // Clone constructor.
    OpalCCollimator(const std::string& name, OpalCCollimator* parent);
    ParticleMatterInteraction* parmatint_m;
};

#endif // OPAL_OpalCCollimator_HH
