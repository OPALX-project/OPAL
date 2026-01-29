//
// Class OpalFlexibleCollimator
//   The Flexible Collimator element.
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
#ifndef OPAL_OpalFlexibleCollimator_HH
#define OPAL_OpalFlexibleCollimator_HH

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;

class OpalFlexibleCollimator: public OpalElement {

public:

    /// The attributes of class OpalFlexibleCollimator.
    enum {
        FNAME = COMMON,  // The horizontal half-size.
        DESC,
        DUMP,
        SIZE
    };

    /// Exemplar constructor.
    OpalFlexibleCollimator();

    virtual ~OpalFlexibleCollimator();

    /// Make clone.
    virtual OpalFlexibleCollimator* clone(const std::string& name);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalFlexibleCollimator(const OpalFlexibleCollimator&);
    void operator=(const OpalFlexibleCollimator&);

    // Clone constructor.
    OpalFlexibleCollimator(const std::string& name, OpalFlexibleCollimator* parent);

    ParticleMatterInteraction* parmatint_m;
};

#endif // OPAL_OpalFlexibleCollimator_HH
