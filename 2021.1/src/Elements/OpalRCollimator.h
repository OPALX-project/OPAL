//
// Class OpalRCollimator
//   The RCOLLIMATOR element.
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
#ifndef OPAL_OpalRCollimator_HH
#define OPAL_OpalRCollimator_HH

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;

class OpalRCollimator: public OpalElement {

public:

    /// The attributes of class OpalRCollimator.
    enum {
        XSIZE = COMMON,  // The horizontal half-size.
        YSIZE,           // The vertical half-size.
        SIZE
    };

    /// Exemplar constructor.
    OpalRCollimator();

    virtual ~OpalRCollimator();

    /// Make clone.
    virtual OpalRCollimator* clone(const std::string& name);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalRCollimator(const OpalRCollimator&);
    void operator=(const OpalRCollimator&);

    // Clone constructor.
    OpalRCollimator(const std::string& name, OpalRCollimator* parent);

    ParticleMatterInteraction* parmatint_m;
};

#endif // OPAL_OpalRCollimator_HH
